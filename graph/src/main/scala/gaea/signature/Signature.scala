package gaea.signature

import gaea.graph._
import gaea.gene.Gene
import gaea.collection.Collection._

import ladder.statistics.Statistics
import org.apache.commons.math3.stat.inference._

import gremlin.scala._
import org.apache.tinkerpop.gremlin.process.traversal.Order
import org.apache.tinkerpop.gremlin.process.traversal.P._

import scalaz._, Scalaz._
import argonaut._, Argonaut._

object Signature {
  val Intercept = Key[Double]("intercept")
  val Expressions = Key[String]("expressions")
  val Coefficient = Key[Double]("coefficient")
  val Coefficients = Key[String]("coefficients")
  val Level = Key[Double]("level")
  val SampleType = Key[String]("sampleType")

  val signatureStep = StepLabel[Vertex]()
  val expressionStep = StepLabel[Vertex]()
  val individualStep = StepLabel[Vertex]()
  val levelStep = StepLabel[Edge]()
  val gidStep = StepLabel[String]()

  val emptyMap = Map[String, Double]()

  val ks = new KolmogorovSmirnovTest()

  lazy val findBackground: GaeaGraph => Map[String, Seq[Double]] = memoize { graph =>
    signatureBackground(graph)
  }

  def dehydrateCoefficients(vertex: Vertex) (key: String): Map[String, Double] = {
    val raw = vertex.property(key).orElse("")
    Parse.parseOption(raw).map(_.as[Map[String, Double]].getOr(emptyMap)).getOrElse(emptyMap)
  }

  def findSignatures(graph: GaeaGraph): List[Tuple2[Vertex, Map[String, Double]]] = {
    val signatureVertexes = graph.typeVertexes("linearSignature")
    signatureVertexes.map((vertex) => (vertex, dehydrateCoefficients(vertex) ("coefficients")))
  }

  def linkSignaturesToGenes(graph: GaeaGraph): List[Tuple2[Vertex, Map[String, Double]]] = {
    val signatures = findSignatures(graph)
    for ((signatureVertex, coefficients) <- signatures) {
      for ((gene, coefficient) <- coefficients) {
        val geneVertex = Gene.findGene(graph) ("gene:" + gene)
        signatureVertex --- ("hasCoefficient", Coefficient -> coefficient) --> geneVertex
      }

      graph.commit()
    }

    signatures
  }

  def signatureLevel
    (genes: Vector[String])
    (coefficients: Vector[Double])
    (intercept: Double)
    (expressions: Map[String, Double])
      : Double = {

    val relevantGenes = selectKeys[String, Double](expressions) (genes) (0.0)
    val (names, levels) = splitMap[String, Double](relevantGenes)

    dotProduct(coefficients) (intercept) (levels)
  }

  def signatureAppliesTo
    (genes: Vector[String])
    (coefficients: Vector[Double])
    (intercept: Double)
    (threshold: Double)
    (expression: Tuple2[Vertex, Map[String, Double]])
      : Boolean = {

    val (vertex, expressions) = expression
    signatureLevel(genes) (coefficients) (intercept) (expressions) > threshold
  }

  def applyExpressionToSignatures
    (graph: GaeaGraph)
    (expressionVertex: Vertex)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : GaeaGraph = {

    val levels = dehydrateCoefficients(expressionVertex) ("expressions")
    val normalized = Statistics.exponentialNormalization(levels)

    for (signature <- signatures) {
      val (signatureVertex, coefficients) = signature
      val intercept = signatureVertex.property(Intercept).orElse(0.0)
      val (genes, values) = splitMap[String, Double](coefficients)
      val level = signatureLevel(genes) (values) (intercept) (normalized)

      signatureVertex --- ("appliesTo", Level -> level) --> expressionVertex
    }

    graph.commit()
    graph
  }

  def applyExpressionsToSignatures
    (graph: GaeaGraph)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : GaeaGraph = {

    val expressionVertexes = graph.typeVertexes("geneExpression")

    for (expressionVertex <- expressionVertexes) {
      applyExpressionToSignatures(graph) (expressionVertex) (signatures)
    }

    graph
  }

  def signatureCorrelation(graph: GaeaGraph) (a: String) (b: String): Tuple3[Vertex, Vertex, Double] = {
    val query = graph.V.has(Gid, within(List(a, b):_*)).as(signatureStep)
      .outE("appliesTo").as(levelStep)
      .inV.as(expressionStep)
      .select((signatureStep, levelStep, expressionStep))
      .map(q => (q._1.property("gid").orElse(""),
        q._1,
        q._2.property("level").orElse(0.0),
        q._3.property("gid").orElse(""))).toSet

    val signatures = query.groupBy(_._1)
    val aNames = signatures(a).map(_._4)
    val bNames = signatures(b).map(_._4)
    val intersect = aNames.intersect(bNames)

    val levels = signatures.map { kv =>
      val (signatureName, tuple) = kv
      (signatureName, tuple.toArray.filter(t => intersect.contains(t._4)).sortBy(_._4))
    }.toMap

    val score = Statistics.pearson(
      breeze.linalg.Vector[Double](levels(a).map(_._3)),
      breeze.linalg.Vector[Double](levels(b).map(_._3)))

    (levels(a).head._2, levels(b).head._2, score)
  }

  def applySignatureCorrelation(graph: GaeaGraph) (a: String) (b:String): Double = {
    val (vertexA, vertexB, score) = signatureCorrelation(graph) (a) (b)
    vertexA <-- ("correlatesTo") --> vertexB
    graph.commit()
    score
  }

  def correlateAllSignatures(graph: GaeaGraph): GaeaGraph = {
    val signatureNames = graph.V.hasLabel("type")
      .has(Gid, "type:linearSignature")
      .out("hasInstance")
      .toSet.map(_.property("gid").orElse(""))

    val pairs = distinctPairs(signatureNames)
    for ((a, b) <- pairs) {
      val score = applySignatureCorrelation(graph) (a) (b)
      println(a.toString + " <--> " + b.toString + ": " + score)
    }

    graph
  }

  def extractLevels(levels: Seq[Tuple2[Vertex, Edge]]): Map[String, Seq[Double]] = {
    levels.groupBy(a => a._1.property("gid").orElse("")).map { s =>
      (s._1, s._2.map(_._2.property("level").orElse(0.0)).sorted)
    }.toMap
  }

  def signatureBackground(graph: GaeaGraph): Map[String, Seq[Double]] = {
    val levelPairs = graph.typeQuery("geneExpression")
      .inE("appliesTo").as(levelStep)
      .outV.as(signatureStep)
      .select((signatureStep, levelStep))
      .toList

    extractLevels(levelPairs)
  }

  // Eventually filter out these variantClassification values: List("5'Flank", "IGR", "Silent", "Intron")`

  def variantLevels(graph: GaeaGraph) (genes: Seq[String]): Map[String, Seq[Double]] = {
    val levelPairs = Gene.synonymsQuery(graph) (genes)
      .in("inGene")
      .out("effectOf")
      .out("tumorSample")
      .in("expressionFor").as(expressionStep)
      .inE("appliesTo").as(levelStep)
      .outV.as(signatureStep)
      .select((signatureStep, levelStep))
      .toSet.toList

    extractLevels(levelPairs)
  }

  def variantSignificance(graph: GaeaGraph) (genes: Seq[String]): Map[String, Double] = {
    val variants = variantLevels(graph) (genes)
    variants.map { variant =>
      val geneLevels = variant._2
      val background = findBackground(graph)
      val back = background(variant._1)
      val backgroundLevels = shear[Double](geneLevels, back)
      val p = ks.kolmogorovSmirnovTest(backgroundLevels.toArray, geneLevels.toArray, true)

      println("background: " + backgroundLevels.size + " - first: " + backgroundLevels.head + " - levels: " + geneLevels.size + " - total: " + back.toSet.size + " - shorn: " + back.toSet.diff(geneLevels.toSet).size)

      (variant._1, p)
    }
  }

  def highestScoringSamples
    (graph: GaeaGraph)
    (signatures: Seq[String])
    (limit: Long)
    (order: Order)
      : Set[Tuple3[Vertex, Vertex, Vertex]] = {

    graph.V.hasLabel("linearSignature")
      .has(Gid, within(signatures:_*)).as(signatureStep)
      .outE("appliesTo").orderBy("level", order).limit(limit)
      .inV.as(expressionStep)
      .out("expressionFor")
      .has(SampleType, "tumor")
      .out("sampleOf").as(individualStep)
      .select((signatureStep, expressionStep, individualStep)).toSet
  }

  def individualScores
    (graph: GaeaGraph)
    (individuals: Seq[String])
    (signatures: Seq[String])
      : Set[Tuple3[Vertex, Vertex, Edge]] = {

    graph.V.hasLabel("individual")
      .has(Gid, within(individuals.toSeq:_*)).as(individualStep)
      .in("sampleOf").has(SampleType, "tumor")
      .in("expressionFor")
      .inE("appliesTo").as(levelStep)
      .outV.has(Gid, within(signatures:_*)).as(signatureStep)
      .select((signatureStep, individualStep, levelStep)).toSet
  }

  // DEPRECATED as inefficient way to do things -----------------------------------------
  def applySignatureToExpressions
    (graph: GaeaGraph)
    (signature: Tuple2[Vertex, Map[String, Double]])
    (expressions: List[Tuple2[Vertex, Map[String, Double]]])
      : GaeaGraph = {

    val (signatureVertex, coefficients) = signature
    val (genes, values) = splitMap[String, Double](coefficients)
    val intercept = signatureVertex.property(Intercept).orElse(0.0)

    val relevant = expressions.filter(signatureAppliesTo(genes) (values) (intercept) (0.5))

    for (expression <- relevant) {
      val (expressionVertex, levels) = expression
      signatureVertex --- ("appliesTo") --> expressionVertex
    }

    graph.commit()
    graph
  }

  def applySignaturesToExpressions
    (graph: GaeaGraph)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : GaeaGraph = {

    val expressionVertexes = graph.typeVertexes("geneExpression")
    val expressions = expressionVertexes.map((vertex) => (vertex, dehydrateCoefficients(vertex) ("expressions")))

    for (signature <- signatures) {
      applySignatureToExpressions(graph) (signature) (expressions)
    }

    graph
  }
}

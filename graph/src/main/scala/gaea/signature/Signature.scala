package gaea.signature

import gaea.feature.Feature
import gaea.titan.Titan
import gaea.math.Stats
import gaea.collection.Collection._

import gremlin.scala._
import com.thinkaurelius.titan.core.TitanGraph
import org.apache.tinkerpop.gremlin.process.traversal.Order
import org.apache.tinkerpop.gremlin.process.traversal.P._

import scalaz._, Scalaz._
import argonaut._, Argonaut._

object Signature {
  val Name = Key[String]("name")
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
  val nameStep = StepLabel[String]()

  val emptyMap = Map[String, Double]()

  def dehydrateCoefficients(vertex: Vertex) (key: String): Map[String, Double] = {
    val raw = vertex.property(key).orElse("")
    Parse.parseOption(raw).map(_.as[Map[String, Double]].getOr(emptyMap)).getOrElse(emptyMap)
  }

  def findSignatures(graph: TitanGraph): List[Tuple2[Vertex, Map[String, Double]]] = {
    val signatureVertexes = Titan.typeVertexes(graph) ("linearSignature")
    signatureVertexes.map((vertex) => (vertex, dehydrateCoefficients(vertex) ("coefficients")))
  }

  def linkSignaturesToFeatures(graph: TitanGraph): List[Tuple2[Vertex, Map[String, Double]]] = {
    val signatures = findSignatures(graph)
    for ((signatureVertex, coefficients) <- signatures) {
      for ((feature, coefficient) <- coefficients) {
        val featureVertex = Feature.findFeature(graph) ("feature:" + feature)
        signatureVertex --- ("hasCoefficient", Coefficient -> coefficient) --> featureVertex
      }

      graph.tx.commit()
    }

    signatures
  }

  def signatureLevel
    (features: Vector[String])
    (coefficients: Vector[Double])
    (intercept: Double)
    (expressions: Map[String, Double])
      : Double = {

    val relevantFeatures = selectKeys[String, Double](expressions) (features) (0.0)
    val (genes, levels) = splitMap[String, Double](relevantFeatures)

    dotProduct(coefficients) (intercept) (levels)
  }

  def signatureAppliesTo
    (features: Vector[String])
    (coefficients: Vector[Double])
    (intercept: Double)
    (threshold: Double)
    (expression: Tuple2[Vertex, Map[String, Double]])
      : Boolean = {

    val (vertex, expressions) = expression
    signatureLevel(features) (coefficients) (intercept) (expressions) > threshold
  }

  def applyExpressionToSignatures
    (graph: TitanGraph)
    (expressionVertex: Vertex)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : TitanGraph = {

    val levels = dehydrateCoefficients(expressionVertex) ("expressions")
    val normalized = Stats.exponentialNormalization(levels)

    for (signature <- signatures) {
      val (signatureVertex, coefficients) = signature
      val intercept = signatureVertex.property(Intercept).orElse(0.0)
      val (features, values) = splitMap[String, Double](coefficients)
      val level = signatureLevel(features) (values) (intercept) (normalized)

      signatureVertex --- ("appliesTo", Level -> level) --> expressionVertex
    }

    graph.tx.commit()
    graph
  }

  def applyExpressionsToSignatures
    (graph: TitanGraph)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : TitanGraph = {

    val expressionVertexes = Titan.typeVertexes(graph) ("geneExpression")

    for (expressionVertex <- expressionVertexes) {
      applyExpressionToSignatures(graph) (expressionVertex) (signatures)
    }

    graph
  }

  def signatureCorrelation(graph: TitanGraph) (a: String) (b: String): Tuple3[Vertex, Vertex, Double] = {
    val query = graph.V.has(Name, within(List(a, b):_*)).as(signatureStep)
      .outE("appliesTo").as(levelStep)
      .inV.as(expressionStep)
      .select((signatureStep, levelStep, expressionStep))
      .map(q => (q._1.property("name").orElse(""),
        q._1,
        q._2.property("level").orElse(0.0),
        q._3.property("name").orElse(""))).toSet

    val signatures = query.groupBy(_._1)
    val aNames = signatures(a).map(_._4)
    val bNames = signatures(b).map(_._4)
    val intersect = aNames.intersect(bNames)

    val levels = signatures.map { kv =>
      val (signatureName, tuple) = kv
      (signatureName, tuple.toArray.filter(t => intersect.contains(t._4)).sortBy(_._4))
    }.toMap

    val score = Stats.pearson(
      breeze.linalg.Vector[Double](levels(a).map(_._3)),
      breeze.linalg.Vector[Double](levels(b).map(_._3)))

    (levels(a).head._2, levels(b).head._2, score)
  }

  def applySignatureCorrelation(graph: TitanGraph) (a: String) (b:String): Double = {
    val (vertexA, vertexB, score) = signatureCorrelation(graph) (a) (b)
    vertexA <-- ("correlatesTo") --> vertexB
    graph.tx.commit()
    score
  }

  def correlateAllSignatures(graph: TitanGraph): TitanGraph = {
    val signatureNames = graph.V.hasLabel("type")
      .has(Name, "type:linearSignature")
      .out("hasInstance")
      .toSet.map(_.property("name").orElse(""))

    val pairs = distinctPairs(signatureNames)
    for ((a, b) <- pairs) {
      val score = applySignatureCorrelation(graph) (a) (b)
      println(a.toString + " <--> " + b.toString + ": " + score)
    }

    graph
  }

  def highestScoringSamples
    (graph: TitanGraph)
    (signatures: Seq[String])
    (limit: Long)
    (order: Order)
      : Set[Tuple3[Vertex, Vertex, Vertex]] = {

    graph.V.hasLabel("linearSignature")
      .has(Name, within(signatures:_*)).as(signatureStep)
      .outE("appliesTo").orderBy("level", order).limit(limit)
      .inV.as(expressionStep)
      .out("expressionFor")
      .has(SampleType, "tumor")
      .out("sampleOf").as(individualStep)
      .select((signatureStep, expressionStep, individualStep)).toSet
  }

  def individualScores
    (graph: TitanGraph)
    (individuals: Seq[String])
    (signatures: Seq[String])
      : Set[Tuple3[Vertex, Vertex, Edge]] = {

    graph.V.hasLabel("individual")
      .has(Name, within(individuals.toSeq:_*)).as(individualStep)
      .in("sampleOf").has(SampleType, "tumor")
      .in("expressionFor")
      .inE("appliesTo").as(levelStep)
      .outV.has(Name, within(signatures:_*)).as(signatureStep)
      .select((signatureStep, individualStep, levelStep)).toSet
  }

  // DEPRECATED as inefficient way to do things -----------------------------------------
  def applySignatureToExpressions
    (graph: TitanGraph)
    (signature: Tuple2[Vertex, Map[String, Double]])
    (expressions: List[Tuple2[Vertex, Map[String, Double]]])
      : TitanGraph = {

    val (signatureVertex, coefficients) = signature
    val (features, values) = splitMap[String, Double](coefficients)
    val intercept = signatureVertex.property(Intercept).orElse(0.0)

    val relevant = expressions.filter(signatureAppliesTo(features) (values) (intercept) (0.5))

    for (expression <- relevant) {
      val (expressionVertex, levels) = expression
      signatureVertex --- ("appliesTo") --> expressionVertex
    }

    graph.tx.commit()
    graph
  }

  def applySignaturesToExpressions
    (graph: TitanGraph)
    (signatures: List[Tuple2[Vertex, Map[String, Double]]])
      : TitanGraph = {

    val expressionVertexes = Titan.typeVertexes(graph) ("geneExpression")
    val expressions = expressionVertexes.map((vertex) => (vertex, dehydrateCoefficients(vertex) ("expressions")))

    for (signature <- signatures) {
      applySignatureToExpressions(graph) (signature) (expressions)
    }

    graph
  }
}

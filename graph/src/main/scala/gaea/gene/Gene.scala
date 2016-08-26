package gaea.gene

import gaea.titan.Titan
import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._
import gaea.collection.Collection._
import org.apache.tinkerpop.gremlin.process.traversal.P._

object Gene {
  val Gid = Key[String]("gid")
  val synonymPrefix = "geneSynonym:"

  val individualStep = StepLabel[Vertex]()
  val variantStep = StepLabel[Vertex]()
  val geneStep = StepLabel[Vertex]()

  def synonymQuery(graph: TitanGraph) (name: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("geneSynonym").has(Gid, synonymPrefix + name).out("synonymFor")
  }

  def synonymsQuery(graph: TitanGraph) (names: Seq[String]): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("geneSynonym").has(Gid, within(names.map(synonymPrefix + _):_*)).out("synonymFor")
  }

  def findSynonymVertex(graph: TitanGraph) (name: String): Option[Vertex] = {
    synonymQuery(graph) (name).headOption
  }

  def findSynonym(graph: TitanGraph) (name: String): Option[String] = {
    val values = findSynonymVertex(graph) (name).map(_.valueMap())
    values.map(vertex => Titan.removePrefix(vertex("gid").asInstanceOf[String]))
  }

  def findGene(graph: TitanGraph) (name: String): Vertex = {
    val inner = Titan.removePrefix(name)
    val synonymName = "geneSynonym:" + inner
    val geneName = "gene:" + inner
    findSynonymVertex(graph) (inner).getOrElse {
      val synonym = graph.V.hasLabel("geneSynonym").has(Gid, synonymName).headOption.getOrElse {
        graph + ("geneSynonym", Gid -> synonymName)
      }

      val gene = graph.V.hasLabel("gene").has(Gid, geneName).headOption.getOrElse {
        graph + ("gene", Gid -> geneName)
      }

      synonym --- ("synonymFor") --> gene
      gene
    }
  }

  def findIndividualsWithVariants(graph: TitanGraph) (gene: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V
      .hasLabel("gene")
      .has(Gid, "gene:" + gene)
      .in("inGene")
      .out("effectOf")
      .out("tumorSample")
      .out("sampleOf")
  }

  def findTumors(graph: TitanGraph) (gene: String): List[String] = {
    findIndividualsWithVariants(graph) (gene)
      .value[String]("submittedTumorSite")
      .toList
  }

  def findTumorCounts(graph: TitanGraph) (gene: String): Map[String, Int] = {
    groupCount[String](findTumors(graph) (gene))
  }

  def findVariantsForIndividuals(graph: TitanGraph) (individuals: Seq[String]) (genes: Seq[String]): Seq[Tuple3[String, String, String]] = {
    val query = graph.V.hasLabel("gene")
      .has(Gid, within(genes.map(synonymPrefix + _):_*)).as(geneStep)
      .in("inGene")
      .out("effectOf").as(variantStep)
      .out("tumorSample")
      .out("sampleOf")
      .has(Gid, within(individuals:_*)).as(individualStep)
      .select((individualStep, variantStep, geneStep))
      .toList

    query.map { q =>
      val (individual, variant, gene) = q
      (individual.property("gid").orElse(""),
        variant.property("variantType").orElse(""),
        gene.property("gid").orElse(""))
    }
  }
}


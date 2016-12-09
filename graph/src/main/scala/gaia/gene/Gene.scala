package gaia.gene

import gaia.graph._
import gremlin.scala._
import gaia.collection.Collection._
import org.apache.tinkerpop.gremlin.process.traversal.P._

object Gene {
  val synonymPrefix = "geneSynonym:"

  val individualStep = StepLabel[Vertex]()
  val variantStep = StepLabel[Vertex]()
  val geneStep = StepLabel[Vertex]()

  def synonymQuery(graph: GaiaGraph) (name: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("geneSynonym").has(Gid, synonymPrefix + name).out("synonymFor")
  }

  def synonymsQuery(graph: GaiaGraph) (names: Seq[String]): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("geneSynonym").has(Gid, within(names.map(synonymPrefix + _):_*)).out("synonymFor")
  }

  def findSynonymVertex(graph: GaiaGraph) (name: String): Option[Vertex] = {
    synonymQuery(graph) (name).headOption
  }

  def findSynonym(graph: GaiaGraph) (name: String): Option[String] = {
    val values = findSynonymVertex(graph) (name).map(_.valueMap())
    values.map(vertex => Gid.removePrefix(vertex("gid").asInstanceOf[String]))
  }

  def findGene(graph: GaiaGraph) (name: String): Vertex = {
    val inner = Gid.removePrefix(name)
    val synonymName = "geneSynonym:" + inner
    val geneName = "gene:" + inner
    findSynonymVertex(graph) (inner).getOrElse {
      val synonym = graph.V.hasLabel("geneSynonym").has(Gid, synonymName).headOption.getOrElse {
        graph.graph + ("geneSynonym", Gid -> synonymName)
      }

      val gene = graph.V.hasLabel("gene").has(Gid, geneName).headOption.getOrElse {
        graph.graph + ("gene", Gid -> geneName)
      }

      synonym --- ("synonymFor") --> gene
      gene
    }
  }

  def findIndividualsWithVariants(graph: GaiaGraph) (gene: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V
      .hasLabel("gene")
      .has(Gid, "gene:" + gene)
      .in("inGene")
      .out("effectOf")
      .out("tumorSample")
      .out("sampleOf")
  }

  def findTumors(graph: GaiaGraph) (gene: String): List[String] = {
    findIndividualsWithVariants(graph) (gene)
      .value[String]("submittedTumorSite")
      .toList
  }

  def findTumorCounts(graph: GaiaGraph) (gene: String): Map[String, Int] = {
    groupCount[String](findTumors(graph) (gene))
  }

  def findVariantsForIndividuals(graph: GaiaGraph) (individuals: Seq[String]) (genes: Seq[String]): Seq[Tuple3[String, String, String]] = {
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


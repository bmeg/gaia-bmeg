package gaea.feature

import gaea.titan.Titan
import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._
import gaea.collection.Collection._
import org.apache.tinkerpop.gremlin.process.traversal.P._

object Feature {
  val Name = Key[String]("name")
  val synonymPrefix = "featureSynonym:"

  val individualStep = StepLabel[Vertex]()
  val variantStep = StepLabel[Vertex]()
  val featureStep = StepLabel[Vertex]()

  def synonymQuery(graph: TitanGraph) (name: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("featureSynonym").has(Name, synonymPrefix + name).out("synonymFor")
  }

  def synonymsQuery(graph: TitanGraph) (names: Seq[String]): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V.hasLabel("featureSynonym").has(Name, within(names.map(synonymPrefix + _):_*)).out("synonymFor")
  }

  def findSynonymVertex(graph: TitanGraph) (name: String): Option[Vertex] = {
    synonymQuery(graph) (name).headOption
  }

  def findSynonym(graph: TitanGraph) (name: String): Option[String] = {
    val values = findSynonymVertex(graph) (name).map(_.valueMap())
    values.map(vertex => Titan.removePrefix(vertex("name").asInstanceOf[String]))
  }

  def findFeature(graph: TitanGraph) (name: String): Vertex = {
    val inner = Titan.removePrefix(name)
    val synonymName = "featureSynonym:" + inner
    val featureName = "feature:" + inner
    findSynonymVertex(graph) (inner).getOrElse {
      val synonym = graph.V.hasLabel("featureSynonym").has(Name, synonymName).headOption.getOrElse {
        graph + ("featureSynonym", Name -> synonymName)
      }

      val feature = graph.V.hasLabel("feature").has(Name, featureName).headOption.getOrElse {
        graph + ("feature", Name -> featureName)
      }

      synonym --- ("synonymFor") --> feature
      feature
    }
  }

  def findIndividualsWithVariants(graph: TitanGraph) (feature: String): GremlinScala[Vertex, shapeless.HNil] = {
    graph.V
      .hasLabel("feature")
      .has(Name, "feature:" + feature)
      .in("inFeature")
      .out("effectOf")
      .out("tumorSample")
      .out("sampleOf")
  }

  def findTumors(graph: TitanGraph) (feature: String): List[String] = {
    findIndividualsWithVariants(graph) (feature)
      .value[String]("submittedTumorSite")
      .toList
  }

  def findTumorCounts(graph: TitanGraph) (feature: String): Map[String, Int] = {
    groupCount[String](findTumors(graph) (feature))
  }

  def findVariantsForIndividuals(graph: TitanGraph) (individuals: Seq[String]) (genes: Seq[String]): Seq[Tuple3[String, String, String]] = {
    val query = graph.V.hasLabel("feature")
      .has(Name, within(genes.map(synonymPrefix + _):_*)).as(featureStep)
      .in("inFeature")
      .out("effectOf").as(variantStep)
      .out("tumorSample")
      .out("sampleOf")
      .has(Name, within(individuals:_*)).as(individualStep)
      .select((individualStep, variantStep, featureStep))
      .toList

    query.map { q =>
      val (individual, variant, feature) = q
      (individual.property("name").orElse(""),
        variant.property("variantType").orElse(""),
        feature.property("name").orElse(""))
    }
  }
}


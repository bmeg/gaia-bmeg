package gaea.facet

import gaea.titan.Titan
import gaea.feature.Feature
import gaea.signature.Signature
import gaea.collection.Collection._

import org.http4s._
import org.http4s.server._
import org.http4s.dsl._

import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._
import org.apache.tinkerpop.gremlin.process.traversal.Order
import org.apache.tinkerpop.gremlin.process.traversal.P._

import com.typesafe.scalalogging._
import _root_.argonaut._, Argonaut._
import org.http4s.argonaut._

import scala.collection.JavaConversions._

object SignatureFacet extends LazyLogging {
  lazy val graph = Titan.connection
  val Name = Key[String]("name")

  val ilog2 = 1.0 / scala.math.log(2)

  def log2(x: Double): Double = {
    scala.math.log(x) * ilog2
  }

  def coefficientsToJson(coefficients: Map[String, Double]) (key: String) (value: String): Json = {
    coefficients.foldLeft(jEmptyArray) { (json, coefficient) =>
      val (feature, level) = coefficient
      val pair = (key, jString(feature)) ->: (value, jNumber(level).getOrElse(jZero)) ->: jEmptyObject
      pair -->>: json
    }
  }

  def propertiesToJson(properties: Map[String, Any]) (key: String) (value: String): Json = {
    properties.foldLeft(jEmptyArray) { (json, property) =>
      val (name, attribute) = property
      val pair = (key, jString(name)) ->: (value, jString(attribute.toString)) ->: jEmptyObject
      pair -->>: json
    }
  }

  def eventMetadata(eventID: String, eventType: String, datatype: String, weights: Map[String, Double]): Json = {
    val weightsJson = coefficientsToJson(weights) ("feature") ("weight")
    ("eventID", jString(eventID)) ->: ("eventType", jString(eventType)) ->: ("datatype", jString(datatype)) ->: ("featureWeights", weightsJson) ->: jEmptyObject
  }

  def signatureToJson(featureNames: List[String]) (vertex: Vertex): Json = {
    val coefficients = Signature.dehydrateCoefficients(vertex) ("coefficients")
    val relevant = selectKeys[String, Double](coefficients) (featureNames) (0.0)
    val score = relevant.values.foldLeft(0.0) ((s, v) => s + Math.abs(v))
    val signatureName = vertex.property(Name).orElse("no name")
    val metadata = eventMetadata(signatureName, "drug sensitivity signature", "NUMERIC", relevant)
    ("score", jNumber(score).getOrElse(jZero)) ->: ("signatureMetadata", metadata) ->: jEmptyObject
  }

  def clinicalEvent(individualVertexes: Seq[Vertex]) (clinicalName: String): Json = {
    val metadata = eventMetadata(clinicalName, "clinical", "STRING", Map[String, Double]())
    val properties = individualVertexes.map(vertex => (vertex.property("name").orElse(""), vertex.property(clinicalName).orElse(""))).toMap
    val json = propertiesToJson(properties) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def expressionEvent(expressions: Seq[Tuple3[String, Vertex, Map[String, Double]]]) (gene: String): Json = {
    val individuals = expressions.map(_._1)
    val coefficients = expressions.map(_._3)
    val metadata = eventMetadata(gene, "mrna_expression", "NUMERIC", Map[String, Double]())
    val expression = coefficients.map(coefficient => log2(coefficient.get(gene).getOrElse(0.0)))
    val properties = individuals.zip(expression).toMap
    val json = coefficientsToJson(properties) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def mutationEvent(mutations: Seq[Tuple3[String, String, String]]) (gene: String): Json = {
    val metadata = eventMetadata(Feature.removePrefix(gene), "mutation call", "STRING", Map[String, Double]())
    val samples = mutations.groupBy(_._1)
    val variants = samples.map { s =>
      val (individual, variants) = s
      (individual, variants.map(_._2).toSet.mkString(","))
    }.toMap

    val json = propertiesToJson(variants) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def levelEvent(levels: Map[String, Double]) (signature: String): Json = {
    val metadata = eventMetadata(signature, "drug sensitivity score", "NUMERIC", Map[String, Double]())
    val json = coefficientsToJson(levels) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def takeHighest(n: Int) (signature: Vertex): List[String] = {
    Signature.dehydrateCoefficients(signature) ("coefficients").toList.sortWith(_._2 > _._2).take(n).map(_._1)
  }

  val service = HttpService {
    case request @ POST -> Root / "gene" =>
      request.as[Json].flatMap { json => 
        val geneNames = json.as[List[String]].getOr(List[String]())
        val featureVertexes = geneNames.map((name: String) => Feature.findSynonymVertex(graph) (name)).flatten
        val featureNames = featureVertexes.map(feature => Titan.removePrefix(feature.property(Name).orElse("")))
        val signatureVertexes = featureVertexes.flatMap(_.in("hasCoefficient").toList).toSet
        val signatureJson = signatureVertexes.map(signatureToJson(featureNames))
        Ok(signatureJson.asJson)
      }

    case request @ POST -> Root / "mutation" =>
      request.as[Json].flatMap { json =>
        val geneNames = json.as[List[String]].getOr(List[String]())
        val featureVertexes = Feature.synonymsQuery(graph) (geneNames).toList
        val featureNames = featureVertexes.map(feature => Feature.removePrefix(feature.property(Name).orElse("")))
        val significance = Signature.variantSignificance(graph) (geneNames).filter(_._2 < 0.05)
        val signatureVertexes = graph.V.has(Name, within(significance.keys.toList:_*)).toList
        val signatureJson = signatureVertexes.map { vertex =>
          significanceToJson(featureNames) (vertex) (significance(vertex.property("name").orElse("")))
        }

        Ok(signatureJson.asJson)
      }

    case request @ POST -> Root / "sample" =>
      request.as[Json].flatMap { json =>
        val metadata = json.as[Map[String, List[Map[String, String]]]].getOr(Map[String, List[Map[String, String]]]())
        val signatureMetadata = metadata("signatureMetadata")
        val expressionMetadata = metadata("expressionMetadata")
        val clinicalEventMetadata = metadata("clinicalEventMetadata")
        val mutationMetadata = metadata("mutationEventMetadata")

        val signatureNames = signatureMetadata.map(_("eventID"))
        val expressionNames = expressionMetadata.map(_("eventID"))
        val clinicalNames = clinicalEventMetadata.map(_("eventID"))
        val mutationNames = mutationMetadata.map(_("eventID"))

        val highestQuery = Signature.highestScoringSamples(graph) (signatureNames) (100) (Order.decr)
        val lowestQuery = Signature.highestScoringSamples(graph) (signatureNames) (100) (Order.incr)
        val query = highestQuery ++ lowestQuery

        val signatureData = query.map(_._1)
        val geneNames = expressionNames ++ signatureData.flatMap(takeHighest(5))

        val individualData = query.map(_._3)
        val individualNames = individualData.map(_.property("name").orElse(""))
        val levelQuery = Signature.individualScores(graph) (individualNames.toList) (signatureNames)

        val expressionData = query.map { q =>
          val (sig, expression, individual) = q
          val coefficients = Signature.dehydrateCoefficients(expression) ("expressions")
          (individual.property("name").orElse(""), expression, coefficients)
        }

        val mutationData = Feature.findVariantsForIndividuals(graph) (individualNames.toList) (mutationNames)
          .groupBy(_._3)

        val levelData = levelQuery.map { q =>
          val (signature, individual, level) = q
          (signature.property("name").orElse(""),
            individual.property("name").orElse(""),
            level.property("level").orElse(0.0))
        }.groupBy(_._1)

        val individualJson = clinicalNames.foldLeft(jEmptyArray) { (json, clinical) =>
          clinicalEvent(individualData.toList) (clinical) -->>: json
        }

        val expressionJson = geneNames.toSet.foldLeft(individualJson) { (json, gene) =>
          expressionEvent(expressionData.toList) (gene) -->>: json
        }

        val mutationJson = mutationData.keys.foldLeft(expressionJson) { (json, gene) =>
          mutationEvent(mutationData(gene).toList) (gene) -->>: json
        }

        val levelJson = levelData.foldLeft(mutationJson) { (json, score) =>
          val (signature, levelTuples) = score
          val levels = levelTuples.map(level => (level._2, level._3)).toMap
          levelEvent(levels) (signature) -->>: json
        }

        Ok(levelJson)
      }
  }
}

package gaea.facet

import gaea.titan.Titan
import gaea.feature.Feature
import gaea.signature.Signature
import gaea.collection.Collection._
import gaea.facet.event.Event

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
  val Gid = Key[String]("gid")

  def takeHighest(n: Int) (signature: Vertex): List[String] = {
    Signature.dehydrateCoefficients(signature) ("coefficients").toList.sortWith(_._2 > _._2).take(n).map(_._1)
  }

  val service = HttpService {
    case request @ POST -> Root / "gene" =>
      request.as[Json].flatMap { json => 
        val geneNames = json.as[List[String]].getOr(List[String]())
        val featureVertexes = geneNames.map((name: String) => Feature.findSynonymVertex(graph) (name)).flatten
        val featureNames = featureVertexes.map(feature => Titan.removePrefix(feature.property(Gid).orElse("")))
        val signatureVertexes = featureVertexes.flatMap(_.in("hasCoefficient").toList).toSet
        val signatureJson = signatureVertexes.map(Event.signatureToJson(featureNames))
        Ok(signatureJson.asJson)
      }

    case request @ POST -> Root / "mutation" =>
      request.as[Json].flatMap { json =>
        val geneNames = json.as[List[String]].getOr(List[String]())
        val featureVertexes = Feature.synonymsQuery(graph) (geneNames).toList
        val featureNames = featureVertexes.map(feature => Titan.removePrefix(feature.property(Gid).orElse("")))
        val significance = Signature.variantSignificance(graph) (geneNames).filter(_._2 < 0.05)
        val signatureVertexes = graph.V.has(Gid, within(significance.keys.toList:_*)).toList
        val signatureJson = signatureVertexes.map { vertex =>
          Event.significanceToJson(featureNames) (vertex) (significance(vertex.property("gid").orElse("")))
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
        val individualNames = individualData.map(_.property("gid").orElse(""))
        val levelQuery = Signature.individualScores(graph) (individualNames.toList) (signatureNames)

        val expressionData = query.map { q =>
          val (sig, expression, individual) = q
          val coefficients = Signature.dehydrateCoefficients(expression) ("expressions")
          (individual.property("gid").orElse(""), expression, coefficients)
        }

        val mutationData = Feature.findVariantsForIndividuals(graph) (individualNames.toList) (mutationNames)
          .groupBy(_._3)

        val levelData = levelQuery.map { q =>
          val (signature, individual, level) = q
          (signature.property("gid").orElse(""),
            individual.property("gid").orElse(""),
            level.property("level").orElse(0.0))
        }.groupBy(_._1)

        val individualJson = clinicalNames.foldLeft(jEmptyArray) { (json, clinical) =>
          Event.clinicalEvent(individualData.toList) (clinical) -->>: json
        }

        val expressionJson = geneNames.toSet.foldLeft(individualJson) { (json, gene) =>
          Event.expressionEvent(expressionData.toList) (gene) -->>: json
        }

        val mutationJson = mutationData.keys.foldLeft(expressionJson) { (json, gene) =>
          Event.mutationEvent(mutationData(gene).toList) (gene) -->>: json
        }

        val levelJson = levelData.foldLeft(mutationJson) { (json, score) =>
          val (signature, levelTuples) = score
          val levels = levelTuples.map(level => (level._2, level._3)).toMap
          Event.levelEvent(levels) (signature) -->>: json
        }

        Ok(levelJson)
      }
  }
}

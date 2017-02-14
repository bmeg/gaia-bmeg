package gaia.facet

import gaia.graph._
import gaia.collection.Collection._
import gaia.facet.event.Event

import org.http4s._
import org.http4s.server._
import org.http4s.dsl._

import gremlin.scala._
import org.apache.tinkerpop.gremlin.process.traversal.Order
import org.apache.tinkerpop.gremlin.process.traversal.P._

import com.typesafe.scalalogging._
import _root_.argonaut._, Argonaut._
import org.http4s.argonaut._
import scala.collection.JavaConversions._

case class IndividualFacet(root: String) extends GaiaFacet with LazyLogging {
  val TumorSite = Key[String]("submittedTumorSite")

  def findIndividualAttributes(graph: GaiaGraph): Set[String] = {
    graph.typeQuery("individual").toList.flatMap(_.valueMap.keys).toSet
  }

  lazy val cacheIndividualAttributes: GaiaGraph => Set[String] = memoize { graph =>
    findIndividualAttributes(graph)
  }

  def individualSurvivalJson(individual: Vertex): Json = {
    val values = individual.valueMap("gid", "vitalStatus", "deathDaysTo", "submittedTumorType")

    val json = ("name", jString(values("gid").asInstanceOf[String])) ->:
      ("status", jString(values("vitalStatus").asInstanceOf[String])) ->:
      ("tumor", jString(values.get("submittedTumorType").getOrElse("unknown").asInstanceOf[String])) ->:
      jEmptyObject

    values.get("deathDaysTo") match {
      case Some(days) => ("days", jNumber(days.asInstanceOf[Long])) ->: json
      case None => json
    }
  }

  def service(graph: GaiaGraph): HttpService = {
    HttpService {
      case GET -> Root / "tumor" / tumorType =>
        Ok(graph.V.has(Gid, "type:individual").out("hasInstance").has(TumorSite, tumorType).value(Gid).toList.asJson)

      case request @ POST -> Root / "survival" =>
        request.as[Json].flatMap { json =>
          val individualNames = json.as[List[String]].getOr(List[String]())
          val individualVertexes = graph.V.hasLabel("individual").has(Gid, within(individualNames:_*)).toList
          val individualJson = individualVertexes.foldLeft(jEmptyArray) {(array, vertex) =>
            individualSurvivalJson(vertex) -->>: array
          }

          Ok(individualJson)
        }

      case GET -> Root / "attributes" =>
        Ok(cacheIndividualAttributes(graph).asJson)

      case request @ POST -> Root / "values" =>
        request.as[Json].flatMap { json =>
          val clinicalNames = json.as[List[String]].getOr(List[String]())
          val individuals = graph.typeQuery("individual").toList // .map(_.valueMap)
          val individualJson = clinicalNames.foldLeft(jEmptyArray) { (json, clinical) =>
            Event.clinicalEvent(individuals) (clinical) -->>: json
          }

          Ok(individualJson)
        }
    }
  }
}

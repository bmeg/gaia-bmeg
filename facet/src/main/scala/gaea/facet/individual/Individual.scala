package gaea.facet

import gaea.titan.Titan
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

object IndividualFacet extends LazyLogging {
  lazy val graph = Titan.connection

  val Name = Key[String]("name")
  val TumorSite = Key[String]("submittedTumorSite")

  def individualSurvivalJson(individual: Vertex): Json = {
    val values = individual.valueMap("name", "vitalStatus", "deathDaysTo", "submittedTumorType")

    val json = ("name", jString(values("name").asInstanceOf[String])) ->:
      ("status", jString(values("vitalStatus").asInstanceOf[String])) ->:
      ("tumor", jString(values.get("submittedTumorType").getOrElse("unknown").asInstanceOf[String])) ->:
      jEmptyObject

    values.get("deathDaysTo") match {
      case Some(days) => ("days", jNumber(days.asInstanceOf[Long])) ->: json
      case None => json
    }
  }

  val service = HttpService {
    case GET -> Root / "tumor" / tumorType =>
      Ok(graph.V.has(Name, "type:individual").out("hasInstance").has(TumorSite, tumorType).value(Name).toList.asJson)

    case request @ POST -> Root / "survival" =>
      request.as[Json].flatMap { json =>
        val individualNames = json.as[List[String]].getOr(List[String]())
        val individualVertexes = graph.V.hasLabel("individual").has(Name, within(individualNames:_*)).toList
        val individualJson = individualVertexes.foldLeft(jEmptyArray) {(array, vertex) =>
          individualSurvivalJson(vertex) -->>: array
        }

        Ok(individualJson)
      }
  }
}

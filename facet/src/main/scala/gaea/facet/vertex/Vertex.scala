package gaea.facet

import gaea.titan.Titan
import gaea.ingest.Ingest
import gaea.collection.Collection._
import gaea.html.VertexHtml

import org.http4s._
import org.http4s.server._
import org.http4s.dsl._
import org.http4s.MediaType._
import org.http4s.headers.{`Content-Type`, `Content-Length`}

import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._
import org.apache.tinkerpop.gremlin.process.traversal.Order
import org.apache.tinkerpop.gremlin.process.traversal.P._

import com.typesafe.scalalogging._
import _root_.argonaut._, Argonaut._
import org.http4s.argonaut._

import scala.collection.JavaConversions._

object VertexFacet extends LazyLogging {
  lazy val graph = Titan.connection
  val Gid = Key[String]("gid")

  def mapToJson(properties: Map[String, Any]) : Json = {
    properties.map( x => {
      ((x._1), (x._2.toString))
    } ).asJson
  }

  def countVertexes(graph: TitanGraph): Map[String, Long] = {
    val counts = graph.V.traversal.label.groupCount.toList.get(0)
    val labels = counts.keySet.toList.asInstanceOf[List[String]]
    labels.foldLeft(Map[String, Long]()) { (countMap, label) =>
      countMap + (label.toString -> counts.get(label))
    }
  }

  lazy val vertexCounts = countVertexes(graph)

  val example = """[{"vertex": "Gene"},
 {"has": {"symbol": ["AHI3", "HOIK4L"]}},
 {"in": "inGene"},
 {"out": "effectOf"},
 {"out": "tumorSample"},
 {"in": "expressionFor"},
 {"as": "expressionStep"},
 {"inE": "appliesTo"},
 {"as": "levelStep"},
 {"outV": ""},
 {"as": "signatureStep"},
 {"select": ["signatureStep", "levelStep", "expressionStep"]}]"""

  val service = HttpService {
    case GET -> Root / "counts" =>
      Ok(vertexCounts.asJson)

    case GET -> Root / "explore" =>
      Ok(VertexHtml.layout(VertexHtml.vertex).toString)
        .withContentType(Some(`Content-Type`(`text/html`)))

    case GET -> Root / "find" / gid =>
      try {
        val vertex = graph.V.has(Gid, gid).head
        val inEdges = groupAs[Edge, String, String](vertex.inE.toList) (_.label) (_.outVertex.value[String]("gid"))
        val outEdges = groupAs[Edge, String, String](vertex.outE.toList) (_.label) (_.inVertex.value[String]("gid"))

        val out = Map[String, Json](
          "type" -> vertex.label.asJson,
          "properties" -> mapToJson(vertex.valueMap),
          "in" -> inEdges.asJson,
          "out" -> outEdges.asJson
        )

        Ok(out.asJson)
      }

      catch {
        case _: Throwable => Ok(Map[String, Json]().asJson)
      }

    case request @ POST -> Root / "query" =>
      request.as[Json].flatMap { query =>
        Ok("yellow".asJson)
      }
  }
}

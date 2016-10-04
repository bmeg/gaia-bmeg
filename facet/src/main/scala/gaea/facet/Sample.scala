package gaea.facet

import gaea.graph._
import gaea.gene.Gene
import gaea.signature.Signature
import gaea.collection.Collection._

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

case class SampleFacet(root: String) extends GaeaFacet with LazyLogging {
  val ignoreMutations = List[String]("5'Flank", "IGR", "Silent", "Intron")

  def service(graph: GaeaGraph) = HttpService {
    case request @ POST -> Root / "variant" / gene =>
      request.as[Json].flatMap { json =>
        Ok("yellow".asJson)
      }
  }
}

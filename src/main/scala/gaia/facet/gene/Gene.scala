package gaia.facet

import gaia.graph._
import gaia.gene.Gene
import gaia.collection.Collection._

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

case class GeneFacet(root: String) extends GaiaFacet with LazyLogging {
  def service(graph: GaiaGraph): HttpService = {
    HttpService {
      case GET -> Root / "hello" / name =>
        Ok(jSingleObject("message", jString(s"Hello, ${name}")))

      case GET -> Root / "synonym" / name =>
        val synonym = Gene.findSynonym(graph) (name).getOrElse {
          "no synonym found"
        }
        Ok(jSingleObject(name, jString(synonym)))

      case GET -> Root / gene / "tumor" / "counts" =>
        Ok(Gene.findTumorCounts(graph) (gene).asJson)
    }
  }
}

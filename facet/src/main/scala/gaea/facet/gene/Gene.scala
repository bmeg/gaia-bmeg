package gaea.facet

import gaea.titan.Titan
import gaea.gene.Gene
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

object GeneFacet extends LazyLogging {
  lazy val graph = Titan.connection

  val service = HttpService {
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

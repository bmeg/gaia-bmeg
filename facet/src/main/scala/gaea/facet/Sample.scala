package gaea.facet

import gaea.titan.Titan
import gaea.gene.Gene
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

object SampleFacet extends LazyLogging {
  lazy val graph = Titan.connection
  val Gid = Key[String]("gid")
  val ignoreMutations = List[String]("5'Flank", "IGR", "Silent", "Intron")

  val service = HttpService {
    case request @ POST -> Root / "variant" / gene =>
      request.as[Json].flatMap { json =>
        Ok("yellow".asJson)
      }
  }
}

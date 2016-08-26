package gaea.titan

import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._
import java.lang.{Long => Llong}

object TitanMigration {
  val indexSpec = Map(
    "positionIndex" -> Map(
      "chromosome" -> classOf[String],
      "strand" -> classOf[String],
      "start" -> classOf[Llong],
      "end" -> classOf[Llong]),

    "idIndex" -> Map("id" -> classOf[String]),
    "gidIndex" -> Map("gid" -> classOf[String]),
    "typeIndex" -> Map("type" -> classOf[String]),
    "tumorIndex" -> Map("submittedTumorSite" -> classOf[String])
  )

  def configuration(): Map[String, String] = {
    Map("storage.cassandra.keyspace" -> "bmeg")
  }

  def migrate(): TitanGraph = {
    val config = Titan.configuration(configuration())
    val graph = Titan.connect(config)
    Titan.makeIndexes(graph) (indexSpec)
    graph
  }
}

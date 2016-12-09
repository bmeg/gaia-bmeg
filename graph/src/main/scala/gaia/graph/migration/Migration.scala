package gaia.graph.migration

import gaia.graph._
import gremlin.scala._
import java.lang.{Long => Llong}

object GaiaBmegMigration {
  val indexSpec = Map(
    "positionIndex" -> Map(
      "chromosome" -> classOf[String],
      "strand" -> classOf[String],
      "start" -> classOf[Llong],
      "end" -> classOf[Llong]),

    "idIndex" -> Map("id" -> classOf[String]),
    "gidIndex" -> Map("gid" -> classOf[String]),
    "symbolIndex" -> Map("symbol" -> classOf[String]),
    "typeIndex" -> Map("type" -> classOf[String]),
    "tumorIndex" -> Map("tumorSite" -> classOf[String])
  )

  def migrate(graph: GaiaGraph): Unit = {
    graph.makeIndexes(indexSpec)
  }
}

package gaia.graph.migration

import gaia.graph._
import gaia.migration._
import gremlin.scala._
import java.lang.{Long => Llong}

object GaiaBmegMigration extends GaiaMigration {
  val indexSpec = Map(
    "positionIndex" -> Map(
      "chromosome" -> classOf[String],
      "strand" -> classOf[String],
      "start" -> classOf[Llong],
      "end" -> classOf[Llong]),

    "symbolIndex" -> Map("symbol" -> classOf[String]),
    "tumorIndex" -> Map("tumorSite" -> classOf[String])
  )

  def signature: String = "gaia-bmeg-migration"

  def migrate(graph: GaiaGraph): Unit = {
    graph.makeIndexes(indexSpec)
  }
}

GaiaMigrations.registerMigrations(List(GaiaBmegMigration))

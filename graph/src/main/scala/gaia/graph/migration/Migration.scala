package gaia.graph.migration

import gaia.graph._
import gremlin.scala._
import java.lang.{Long => Llong}

object GaiaBmegMigration extends GaiaMigration {
  GaiaMigrations.registerMigrations(List(GaiaBmegMigration))

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

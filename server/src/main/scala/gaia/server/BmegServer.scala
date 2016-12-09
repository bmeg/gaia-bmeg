package gaia.server

import gaia.config._
import gaia.graph._
import gaia.facet._

object BmegServer extends App {
  val config = GaiaConfig.readConfig("resources/config/gaia.yaml")
  val graph = config.connectToGraph(config.graph)

  if (graph.isSuccess) {
    GaiaServer.start(config.server) (graph.get) (BmegFacets.facets)
  } else {
    println("failed to connect to graph: " + config.graph.toString)
  }
}

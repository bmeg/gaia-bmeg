package gaea.server

import gaea.config._
import gaea.graph._
import gaea.facet._

object BmegServer extends App {
  val config = GaeaConfig.readConfig("resources/config/gaea.yaml")
  val graph = config.connectToGraph(config.graph)

  if (graph.isSuccess) {
    GaeaServer.start(config.server) (graph.get) (BmegFacets.facets)
  } else {
    println("failed to connect to graph: " + config.graph.toString)
  }
}

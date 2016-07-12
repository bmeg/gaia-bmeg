package gaea.server

import gaea.facet.BmegFacets

object BmegServer extends App {
  GaeaServer.start(BmegFacets.facets)
}

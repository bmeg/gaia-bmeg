package gaea.facet

import org.http4s._

object BmegFacets {
  val facets = List[GaeaFacet](
    new IndividualFacet("/individual/"),
    new GeneFacet("/gene/"),
    new SignatureFacet("/signature/")
  )
}

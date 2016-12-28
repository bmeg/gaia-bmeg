package gaia.facet

import org.http4s._

object BmegFacets {
  val facets = List[GaiaFacet](
    new IndividualFacet("/individual/"),
    new GeneFacet("/gene/"),
    new SignatureFacet("/signature/")
  )
}

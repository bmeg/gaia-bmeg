package gaea.facet

import org.http4s._

object BmegFacets {
  val facets = List[Tuple2[String, HttpService]](
    ("/individual/", IndividualFacet.service),
    ("/feature/", FeatureFacet.service),
    ("/signature/", SignatureFacet.service)
  )
}

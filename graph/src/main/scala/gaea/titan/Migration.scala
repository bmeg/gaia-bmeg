package gaea.titan

import gaea.feature.Hugo

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

    "nameIndex" -> Map("name" -> classOf[String]),
    "genderIndex" -> Map("gender" -> classOf[String]),
    "cancerIndex" -> Map("submittedTumorSite" -> classOf[String])
  )

  def migrate(): TitanGraph = {
    val config = Titan.configuration(Map[String, String]())
    val graph = Titan.connect(config)
    Titan.makeIndexes(graph) (indexSpec)
    Hugo.hugoMigration(graph) ("resources/hugo-names")
    graph
  }

  def signatureMigration(graph: TitanGraph): Vertex = {
    val sigNames = List("linearSignature:17-AAG_median", "linearSignature:AEW541_median", "linearSignature:AZD6244_lowerquartile", "linearSignature:Irinotecan_median", "linearSignature:Paclitaxel_median", "linearSignature:Panobinostat_median", "linearSignature:PD-0325901_lowerquartile", "linearSignature:RAF265_median", "linearSignature:TAE684_median", "linearSignature:TKI258_median", "linearSignature:Topotecan_median")

    val sig = graph + ("type", Key[String]("name") -> "type:linearSignature")
    val sigV = sigNames.map(graph.V.hasLabel("linearSignature").has(Key[String]("name"), _).head)
    for (sigVertex <- sigV) {sig --- ("hasInstance") --> sigVertex}
    graph.tx.commit()
    sig
  }
}

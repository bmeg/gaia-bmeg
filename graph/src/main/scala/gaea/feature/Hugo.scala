package gaea.feature

import gaea.titan.Titan

import scala.io.Source
import com.thinkaurelius.titan.core.TitanGraph
import gremlin.scala._

object Hugo {
  val HugoId = Key[String]("hugoId")
  val Name = Key[String]("name")
  val Description = Key[String]("description")
  val Chromosome = Key[String]("chromosome")
  val Accession = Key[String]("accession")
  val Refseq = Key[String]("refseq")

  def readHugo(filename: String): List[Array[String]] = {
    val hugoLines = Source.fromFile(filename).getLines
    val header = hugoLines.next
    val allHugos = hugoLines.map(_.split("\t")).toList
    allHugos.filter(_(3) == "Approved")
  }

  def featureConvoy(graph: TitanGraph) (featureType: Vertex) (synonymType: Vertex) (hugo: Array[String]): Vertex = {
    val chromosome = if(hugo.length > 6) hugo(6) else ""
    val accession = if(hugo.length > 7) hugo(7) else ""
    val refseq = if(hugo.length > 8) hugo(8) else ""
    val feature = graph + ("feature",
      Name -> ("feature:" + hugo(1)),
      HugoId -> hugo(0).replaceFirst("HGNC:", ""),
      Description -> hugo(2),
      Chromosome -> chromosome,
      Accession -> accession,
      Refseq -> refseq)
    featureType --- ("hasInstance") --> feature

    val otherSynonyms = if(hugo.length > 5 && hugo(5) != "") hugo(5).split(", ") else Array[String]()
    val synonyms = otherSynonyms :+ hugo(1)
    for (synonym <- synonyms) {
      val synonymVertex = graph + ("featureSynonym", Name -> ("featureSynonym:" + synonym))
      synonymType --- ("hasInstance") --> synonymVertex
      synonymVertex --- ("synonymFor") --> feature
    }

    feature
  }

  def hugoConvoy(graph: TitanGraph) (hugos: List[Array[String]]): Integer = {
    val featureType = graph + ("type", Name -> "type:feature")
    val synonymType = graph + ("type", Name -> "type:featureSynonym")

    for (hugo <- hugos) {
      print(".")
      featureConvoy(graph) (featureType) (synonymType) (hugo)
    }

    hugos.length
  }

  def hugoMigration(graph: TitanGraph) (hugoFile: String): Integer = {
    val hugos = readHugo(hugoFile)
    val count = hugoConvoy(graph) (hugos)
    graph.tx.commit()
    count
  }
}


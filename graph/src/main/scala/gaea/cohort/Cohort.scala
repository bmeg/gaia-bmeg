package gaea.cohort

import gaea.graph._
import gremlin.scala._
import gaea.frame.Frame
import gaea.collection.Collection._
import org.apache.tinkerpop.gremlin.process.traversal.P._

object Cohort {
  val sampleStep = StepLabel[Vertex]()
  val geneStep = StepLabel[String]()
  val biosamplePrefix = "biosample:".r

  class VariantBuilder(header: List[String], data: List[List[String]]) {
    def addSample(sample: String) (variants: List[String]): VariantBuilder = {
      val variantSet = variants.toSet
      val newKeys = variantSet.diff(header.toSet).toList
      val newHeader = header ++ newKeys
      val row = sample :: newHeader.map(gene => if (variantSet(gene)) "1" else "0")
      val newData = row :: data
      new VariantBuilder(newHeader, newData)
    }

    def finish(): Seq[Seq[String]] = {
      val lines = ("gene" :: header) :: data
      val headerSize = lines.head.size
      val output = lines.map { line =>
        line ++ List.fill(headerSize - line.size) ("0")
      }

      transpose(output)
    }
  }

  def transpose[A](matrix: Seq[Seq[A]]): Seq[Seq[A]] = {
    val base = List.fill(matrix.head.size) (Vector[A]())
    matrix.foldLeft(base) { (t, row) =>
      t.zip(row).map { case (column, cell) =>
        column :+ cell
      }
    }
  }

  def cohortVariants(graph: GaeaGraph) (cohort: String): Seq[Seq[String]] = {
    val variants = graph.V
      .hasLabel("cohort")
      .has(Gid, "cohort:" + cohort)
      .out("hasMember").as(sampleStep)
      .in("tumorSample")
      .in("effectOf")
      .out("inGene")
      .value[String]("gid").as(geneStep)
      // .value[String]("symbol").as(geneStep)
      .select((sampleStep, geneStep)).toList

    // val samples = groupAs(variants) (_._1.property("gid").orElse("")) (_._2.substring(5))
    val samples = groupAs(variants) (_._1.property("gid").orElse("")) (_._2.substring(5))
    val base = new VariantBuilder(List[String](), List[List[String]]())
    samples.foldLeft(base) { case (builder, (sample, genes)) =>
      builder.addSample(sample) (genes)
    }.finish()
  }

  def cohortVariantsTSV(graph: GaeaGraph) (cohort: String): String = {
    val variants = cohortVariants(graph) (cohort)
    val frame = Frame.renderTSV(variants)
    biosamplePrefix.replaceAllIn(frame, "")
  }
}

package gaea.facet.event

import gaea.graph._
import gaea.signature.Signature
import gaea.collection.Collection._

import _root_.argonaut._, Argonaut._
import org.http4s.argonaut._

import gremlin.scala._

import scala.collection.JavaConversions._

object Event {
  val ilog2 = 1.0 / scala.math.log(2)

  def log2(x: Double): Double = {
    scala.math.log(x) * ilog2
  }

  def coefficientsToJson(coefficients: Map[String, Double]) (key: String) (value: String): Json = {
    coefficients.foldLeft(jEmptyArray) { (json, coefficient) =>
      val (gene, level) = coefficient
      val pair = (key, jString(gene)) ->: (value, jNumber(level).getOrElse(jZero)) ->: jEmptyObject
      pair -->>: json
    }
  }

  def propertiesToJson(properties: Map[String, Any]) (key: String) (value: String): Json = {
    properties.foldLeft(jEmptyArray) { (json, property) =>
      val (gid, attribute) = property
      val pair = (key, jString(gid)) ->: (value, jString(attribute.toString)) ->: jEmptyObject
      pair -->>: json
    }
  }

  def eventMetadata(eventID: String, eventType: String, datatype: String, weights: Map[String, Double]): Json = {
    val weightsJson = coefficientsToJson(weights) ("gene") ("weight")
    ("eventID", jString(eventID)) ->: ("eventType", jString(eventType)) ->: ("datatype", jString(datatype)) ->: ("geneWeights", weightsJson) ->: jEmptyObject
  }

  def signatureToJson(geneNames: List[String]) (vertex: Vertex): Json = {
    val coefficients = Signature.dehydrateCoefficients(vertex) ("coefficients")
    val relevant = selectKeys[String, Double](coefficients) (geneNames) (0.0)
    val score = relevant.values.foldLeft(0.0) ((s, v) => s + Math.abs(v))
    val signatureName = vertex.property(Gid).orElse("no name")
    val metadata = eventMetadata(signatureName, "drug sensitivity signature", "NUMERIC", relevant)
    ("score", jNumber(score).getOrElse(jZero)) ->: ("signatureMetadata", metadata) ->: jEmptyObject
  }

  def significanceToJson(geneNames: List[String]) (vertex: Vertex) (significance: Double): Json = {
    val coefficients = Signature.dehydrateCoefficients(vertex) ("coefficients")
    val relevant = selectKeys[String, Double](coefficients) (geneNames) (0.0)
    val signatureName = vertex.property(Gid).orElse("no name")
    val metadata = eventMetadata(signatureName, "significance to mutation", "NUMERIC", relevant)
    ("significance", jNumber(significance).getOrElse(jZero)) ->: ("signatureMetadata", metadata) ->: jEmptyObject
  }

  def clinicalEvent(individualVertexes: Seq[Vertex]) (clinicalName: String): Json = {
    val metadata = eventMetadata(clinicalName, "clinical", "STRING", Map[String, Double]())
    val properties = individualVertexes.map(vertex => (vertex.property("gid").orElse(""), vertex.property(clinicalName).orElse(""))).toMap
    val json = propertiesToJson(properties) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def expressionEvent(expressions: Seq[Tuple3[String, Vertex, Map[String, Double]]]) (gene: String): Json = {
    val individuals = expressions.map(_._1)
    val coefficients = expressions.map(_._3)
    val metadata = eventMetadata(gene, "mrna_expression", "NUMERIC", Map[String, Double]())
    val expression = coefficients.map(coefficient => log2(coefficient.get(gene).getOrElse(0.0)))
    val properties = individuals.zip(expression).toMap
    val json = coefficientsToJson(properties) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def mutationEvent(mutations: Seq[Tuple3[String, String, String]]) (gene: String): Json = {
    val metadata = eventMetadata(Gid.removePrefix(gene), "mutation call", "STRING", Map[String, Double]())
    val samples = mutations.groupBy(_._1)
    val variants = samples.map { s =>
      val (individual, variants) = s
      (individual, variants.map(_._2).toSet.mkString(","))
    }.toMap

    val json = propertiesToJson(variants) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }

  def levelEvent(levels: Map[String, Double]) (signature: String): Json = {
    val metadata = eventMetadata(signature, "drug sensitivity score", "NUMERIC", Map[String, Double]())
    val json = coefficientsToJson(levels) ("sampleID") ("value")
    ("metadata", metadata) ->: ("sampleData", json) ->: jEmptyObject
  }
}

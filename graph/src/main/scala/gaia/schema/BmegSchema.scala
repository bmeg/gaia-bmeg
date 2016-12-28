package gaia.schema

import org.json4s._
import org.json4s.jackson.JsonMethods._
import org.json4s.jackson.Serialization.{write}

// sealed trait PropertyTranslation[A, Z] {
//   def translate(property: A): Map[String, Z]
//   def key: String
// }

// case class StringTranslation(key: String) extends PropertyTranslation[Any, String] {
//   def translate(property: Any): Map[String, String] = {
//     Map(key -> property.toString)
//   }
// }

// case class SerializeTranslation[P](key: String) extends PropertyTranslation[Map[String, P], String] {
//   def translate(property: Map[String, P]): Map[String, String] = {
//     val serialized = 
//     Map(key -> )
//   }
// }

object BmegSchema {
  val gene = GaiaVertex("gene", "type", Map("symbol" -> "string"))
  val geneSynonym = GaiaVertex("geneSynonym", "type", Map("symbol" -> "string"))
  val geneSynonym->gene = GaiaEdge("geneSynonym", "synonymFor", "gene")

  val vertexes = List(gene, geneSynonym)
  val edges = List(geneSynonym->gene)

  val schema = GaiaSchema.assemble(vertexes, edges)
}

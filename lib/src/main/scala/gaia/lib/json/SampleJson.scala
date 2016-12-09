package gaia.lib.json

import scalaz._, Scalaz._
import argonaut._, Argonaut._
import java.io.File
import scala.util.matching.Regex
import scala.io.Source

object SampleJson {
  def listDir(dir: String): List[File] = {
    val directory = new File(dir)
    if (directory.exists && directory.isDirectory) {
      directory.listFiles.toList
    } else {
      List[File]()
    }
  }

  def listFiles(dir: String): List[File] = {
    listDir(dir).filter(_.isFile)
  }

  def matchFiles(dir: String) (regex: Regex): List[File] = {
    listFiles(dir).filter((file: File) => !regex.findAllIn(file.getName).isEmpty)
  }

  def parseJsonLines(file: File): List[Json] = {
    val lines = Source.fromFile(file).getLines
    lines.map(Parse.parseOption(_)).toList.flatten
  }

  def parseJsonFiles(files: List[File]): List[Json] = {
    files.flatMap(parseJsonLines(_))
  }

  def parseJsonMatch(dir: String) (regex: Regex): List[Json] = {
    val files = matchFiles(dir) (regex)
    parseJsonFiles(files)
  }

  def jsonFields(json: List[Json]) (field: String): List[String] = {
    json.map(_.field(field)).flatten.map(_.toString)
  }
}

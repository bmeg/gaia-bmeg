package gaia.server

import gaia.config._
import gaia.graph._
import gaia.facet._

object BmegServer extends App {
  val parser = new gaia.command.GaiaServerCommand("0.0.7")
  parser.execute(args)
}

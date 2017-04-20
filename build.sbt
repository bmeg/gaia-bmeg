organization := "io.bmeg"
name := "bmeg-graph"
version := "0.0.1-SNAPSHOT"

scalaVersion := "2.11.8"
resolvers += "Local Maven Repository" at "file://"+Path.userHome.absolutePath+"/.m2/repository"

resolvers ++= Seq(
  "Akka Repository" at "http://repo.akka.io/releases/",
  "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots",
  "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases",
  "Twitter Maven Repo" at "https://maven.twttr.com",
  "GAEA Depends Repo" at "https://github.com/bmeg/gaia-depends/raw/master/"
)

libraryDependencies ++= Seq(
  "org.apache.commons" %  "commons-math3"   % "3.6.1",
  "io.bmeg"            %% "gaia-server"     % "0.0.9-SNAPSHOT",
  "io.bmeg"            %% "ladder"          % "0.0.1-SNAPSHOT"
)

libraryDependencies += "org.slf4j" % "slf4j-simple" % "1.6.4"

publishTo := {
  val nexus = "https://oss.sonatype.org/"
  if (isSnapshot.value)
    Some("snapshots" at nexus + "content/repositories/snapshots")
  else
    Some("releases"  at nexus + "content/repositories/releases")
}

credentials += Credentials(Path.userHome / ".ivy2" / ".credentials")

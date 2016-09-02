# gaea-bmeg

BioMedical Evidence Graph for Gaea

## requirements

You must have installed:

* java 8
* cassandra
* zookeeper
* kafka
* sbt

## setup

Start cassandra, zookeeper and kafka. Note the addresses for cassandra and kafka as we will need these later.

Clone this repo, then `cd` into it.

Type the command:

    sbt server/console

We will need to create some indexes... do that now:

```
import gaea.titan.Titan
val graph = Titan.defaultGraph()
Titan.makeIndex(graph) ("gidIndex") (Map("gid" -> classOf[String]))
Titan.makeIndex(graph) ("symbolIndex") (Map("symbol" -> classOf[String]))
```

Then, let's ingest some vertexes. For this example, we will be pulling the pancreatic data from TCGA.

```
val ingestor = gaea.kafka.Ingestor("10.96.11.91:9092", "gaea-bmeg-test", List("gaea-hugo", "gaea-paad"))
ingestor.ingest()
```

Let this run until it is done. You are all set!
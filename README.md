# gaia-bmeg

BioMedical Evidence Graph using the Gaia platform

![BMEG](https://github.com/bmeg/gaia-bmeg/blob/master/resources/public/static/img/schema.png)

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

First we initialize the project. Type the command:

    sbt "server/run init"

This creates some indexes and prepares the database. Next, let's ingest some things:

```
val ingestor = new gaia.kafka.Ingestor("10.96.11.89:9092", "pick-a-unique-group-id", List("hugo.Gene", "ccle.Biosample"))
ingestor.ingest()
```

Let this run until it is done (you will have to kill the process as it will wait forever for additional messages). You are all set!

## additional topics

There are several topics each containing their own data set which can be loaded into Gaia.

* hugo.Gene
* hugo.GeneSynonym
* hugo.GeneFamily
* hugo.GeneDatabase
* ccle.Biosample
* ccle.Cohort
* ccle.GeneExpression
* ccle.ResponseCurve
* ccle.CallSet
* ccle.Variant
* ccle.VariantAnnotation
* ccle.CNACallSet
* ccle.CNASegment
* ctdd.Compound
* ctdd.ResponseCurve
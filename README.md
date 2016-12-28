# gaia-bmeg

BioMedical Evidence Graph for Gaia

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
gaia.titan.TitanMigration.migrate()
```

Then, let's ingest some vertexes. For this example, we will be pulling the pancreatic data from TCGA.

```
val ingestor = new gaia.kafka.Ingestor("10.96.11.91:9092", "pick-a-unique-group-id", List("gaia-hugo", "gaia-paad"))
ingestor.ingest()
```

Let this run until it is done (you will have to kill the process as it will wait forever for additional messages). You are all set!

## additional topics

There are several topics each containing their own data set which can be loaded into Gaia.

* gaia-hugo : All the gene names and their synonyms
* gaia-paad : Just the pancreatic dataset from TCGA (useful for building a minimal working DB)
* gaia-prad : Just the prostate dataset from TCGA (again, a good subset for playing around)
* ccle-phenotypes : All the variant and drug response data from CCLE
* tcga-expression : All the gene expression data from TCGA
* tcga-clinical : Clinical information about the patients in the TCGA
* tcga-variants : All the variant calls and associated samples from TCGA

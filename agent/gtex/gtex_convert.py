#!/usr/bin/env python

# http://www.gtexportal.org/static/datasets/biobank/downloads/biobank_collection_20170329_093753.txt
# http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz

from ga4gh.schemas.ga4gh import bio_metadata_pb2
from bmeg import matrix_pb2
from google.protobuf.json_format import MessageToDict

import csv
import json
import pandas
import argparse

INDIVIDUAL_HEADERS = [
    "hasBrainTissue",
    "pathologyNotes",
    "rin",
    "ageBracket",
    "hardyScale"
]

SAMPLE_HEADERS = [
    "amount",
    "autolysisScore",
    "concentration",
    "expression",
    "genotype",
    "hasGTExImage",
    "hasBRDImage",
    "mass",
    "materialType",
    "originalMaterialType",
    "tissueSite",
    "tissueSiteDetail",
    "tissueSiteDetailId",
    "volume"
]


def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a


def emit_json(prefix, out_handles, message):
    name = message.DESCRIPTOR.full_name
    if name not in out_handles:
        out_handles[name] = open(prefix + name, "w")
    out_handles[name].write(json.dumps(MessageToDict(message)))
    out_handles[name].write("\n")

def parse_bio(bio, out):
    out_handles = {}
    with open(bio) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        
        for row in reader:
            ind = bio_metadata_pb2.Individual()
            ind.id = row['subjectId']
            ind.sex.term = row['sex']
            #out.individual_age_at_collection.age = row['ageBracket']
            for i in INDIVIDUAL_HEADERS:
                proto_list_append( ind.info[i], row[i] )
            emit_json(out, out_handles, ind )

            sam = bio_metadata_pb2.Biosample()
            sam.id = row['sampleId']
            sam.dataset_id = "gtex"
            sam.individual_id = row['subjectId']
            for i in SAMPLE_HEADERS:
                proto_list_append( sam.info[i], row[i] )
            emit_json(out, out_handles, sam)
    for handle in out_handles.values():
        handle.close()
                


def parse_expression(path, out):
    df = pandas.read_csv(path, compression="gzip", sep="\t", index_col=0, skiprows=2)
    
    df = df.groupby("Description").mean().transpose()

    out_handles = {}

    for row in df.iterrows():
        gex = matrix_pb2.GeneExpression()
        gex.id = "GeneExpression:%s" % (row[0])
        gex.source = "gtex"
        gex.biosample_id = "biosample:%s" % (row[0])
        gex.scale = matrix_pb2.RPKM
        for k, v in row[1].iteritems():
            gex.expressions[k] = v
        emit_json(out, out_handles, gex)
    for handle in out_handles.values():
        handle.close()


if __name__ == "__main__":
    # Construct the parser
    parser = argparse.ArgumentParser()

    # Now add all the options to it
    parser.add_argument('--bio', type=str, help='')
    parser.add_argument('--expression', type=str, help='')
    parser.add_argument('--out', default="gtex.")
    args = parser.parse_args()

    if args.bio:
        parse_bio(args.bio, args.out)
    if args.expression:
        parse_expression(args.expression, args.out)
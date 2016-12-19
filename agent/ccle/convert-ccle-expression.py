#!/usr/bin/env python


"""
curl -o CCLE_Expression_2012-09-29.res "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_21/CCLE_Expression_2012-09-29.res?downloadff=true&fileId=6760"
"""

import re
import sys
import csv
import json
from bmeg import phenotype_pb2, sample_pb2, genome_pb2, variant_pb2

from google.protobuf import json_format


def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['type'] = message.DESCRIPTOR.name
    return json.dumps(msg)

with open(sys.argv[1]) as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = reader.next()
    reader.next()
    #reader.next()
    
    values = {}
    for i in header[2:]:
        if len(header):
            values[i] = []

    gene_order = []
    for line in reader:
        if len(line[0]):
            for a in zip(header, line)[2:]:
                if len(a[0]):
                    values[a[0]].append(float(a[1]))
            gene_order.append(line[0])

    with open("expression_ccle.json", "w") as handle:
        for k, v in values.items():
            ge = sample_pb2.GeneExpression()
            ge.gid = "ccle_expression:%s" % (k)
            vals = {}
            counts = {}
            for g, val in zip(gene_order, v):
                if g not in vals:
                    vals[g] = []
                    counts[g] = 0.0
                vals[g].append(val)
                counts[g] += 1.0
            for g in vals:
                ge.expressions[g] = sum(vals[g]) / counts[g]
            handle.write("%s\n" % (message_to_json(ge)))
#!/usr/bin/env python


"""
curl -o CCLE_Expression_2012-09-29.res "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_21/CCLE_Expression_2012-09-29.res?downloadff=true&fileId=6760"
"""


import sys
import csv
import json


with open(sys.argv[1]) as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = reader.next()
    reader.next()
    reader.next()
    
    values = {}
    for i in header[2:]:
        if len(header):
            values[i] = []

    gene_order = []
    for line in reader:
        for a in zip(header, line)[2:]:
            if len(a[0]):
                values[a[0]].append(float(a[1]))
        gene_order.append(line[0])

    with open("expression_ccle.json", "w") as handle:
        for k, v in values.items():
            a = dict(zip(gene_order, v))
            out = {
                "gid" : "ccle_expression:%s" % (k),
                "type" : "GeneExpression",
                "expressions" : a
            }
            handle.write("%s\n" % (json.dumps(out)))
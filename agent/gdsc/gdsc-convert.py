#!/usr/bin/env python

import sys
import pandas
import math
import json
from bmeg import phenotype_pb2
from google.protobuf import json_format
from ga4gh import bio_metadata_pb2

def emit(message):
    msg = json_format.MessageToDict(message)
    msg["#label"] = message.DESCRIPTOR.full_name
    print json.dumps(msg)

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a

def gdsc_ic50_row(row, compound_table, sample_table, emit):
    sample_name = sample_table[row["COSMIC_ID"]]
    compound_name = compound_table[ int(row["DRUG_ID"]) ]
    
    gid = "responseCurve:%s:%s" % (sample_name, compound_name)

    response = phenotype_pb2.ResponseCurve()
    response.gid = gid
    response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
    response.compound = compound_name
    response.sample = sample_name

    s = response.summary.add()
    s.type = phenotype_pb2.ResponseSummary.IC50
    s.value = row['LN_IC50']
    s.unit = "uM"

    s = response.summary.add()
    s.type = phenotype_pb2.ResponseSummary.AUC
    s.value = row['AUC']
    s.unit = "uM"

    emit(response)

def gdsc_cell_info(row, emit):
    sample = bio_metadata_pb2.Biosample()
    sample.id = "biosample:GDSC:%s" % row["Sample Name"]
    sample.dataset_id = "GDSC"

    dis = row['GDSC\nTissue\ndescriptor 2']
    if not isinstance(dis, float):
        sample.disease.term = dis.lower().replace('_', ' ')

    proto_list_append(sample.attributes.attr['sampleType'], "cellline")
    label = row['Cancer Type\n(matching TCGA label)']
    if not isinstance(label, float) and len(label):
        proto_list_append(sample.attributes.attr['source'], label)
    emit(sample)

conv_file = sys.argv[1]
cell_info_file = sys.argv[2]
compound_info_file = sys.argv[3]
fitted_file = sys.argv[4]


cl_info = pandas.read_excel(conv_file, index_col=0)
sample_table = {}
for row in cl_info.iterrows():
    sample_table[row[0]] = "biosample:CCLE:%s" % (row[1]['CCLE name'])

cl_info = pandas.read_excel(cell_info_file, index_col=1)
for row in cl_info.iterrows():
    if row[0] not in sample_table:
        sample_table[row[0]] = "biosample:GDSC:%s" % (row[1]['Sample Name'])
        gdsc_cell_info(row[1], emit)
        

    
compound_table = {}
"""
comp_info = pandas.read_excel(conv_file, sheetname=1, index_col=0)
for row in comp_info.iterrows():
    compound_table[int(row[0])] = row[1]['GDSC name']
"""
comp_info = pandas.read_excel(compound_info_file, index_col=0)
for row in comp_info.iterrows():
    compound_table[int(row[0])] = "compound:%s" % (row[1]['Drug Name'])


    
    
fitted = pandas.read_excel(fitted_file)
for r in fitted.iterrows():
    cosmic_id = int(r[1]["COSMIC_ID"])
    if cosmic_id in cl_info.index:
        gdsc_ic50_row( r[1], compound_table, sample_table, emit )
        
        


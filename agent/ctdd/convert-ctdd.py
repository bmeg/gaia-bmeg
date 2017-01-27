#! /usr/bin/python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

This program converts CTDD drug response information into
protobuf data based on the BMEG sample.proto schema.

Source: ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip

The four files of interest (for this converter) are:
1) v20.data.curves_post_qc.txt
2) v20.meta.per_compound.txt
3) v20.meta.per_cell_line.txt
4) v20.meta.per_experiment.txt
'''

from bmeg import phenotype_pb2
from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re
import pandas

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--response', type=str, help='Path to the drug response experiment data you want to import')
    parser.add_argument('--metadrug', type=str, help='Path to the drug meta data you want to import')
    parser.add_argument('--metacellline', type=str, help='Path to the cell line meta data you want to import')
    parser.add_argument('--metaexperiment', type=str, help='Path to the experiment meta data you want to import')
    parser.add_argument('--data', type=str, help='Path to the experiment data you want to import')
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

########################################

def find_biosample(state, source, barcode, sample_type):
    sample_name = 'biosample:CCLE:' + barcode
    biosample = state['Biosample'].get(sample_name)
    if biosample is None:
        biosample = schema.Biosample()
        biosample.name = sample_name
        biosample.dataset_id = "CCLE"
        biosample.source = source
        biosample.barcode = barcode
        biosample.sampleType = sample_type
        state['Biosample'][sample_name] = biosample

    return biosample

def append_unique(l, i):
    if not i in l:
        l.append(i)

def process_drugs(emit, input): #row is a namedtuple
    drugs = set()
    
    for row in input.itertuples():
        # create drug message for CTDD compound
        drug_name = "drug:" + row.cpd_name # in the future might want to find canonical drug name via external resources
        if drug_name not in drugs:
            drug = phenotype_pb2.Compound()
            drug.id = drug_name
            drug.smiles = row.cpd_smiles
            drug.synonyms.append("broad.org/cpd/" + row.broad_cpd_id)
            #drug.synonyms.append("drug:" + row.cpd_smiles)
            emit(drug)
            drugs.add(drug_name)

def process_response(emit, input, data):
    
    gid_set = set()
    for row in input.itertuples():
        gid = "lincs.org/%s/%s" % (row.ccl_name, row.cpd_name)
        if gid not in gid_set:
            response = phenotype_pb2.ResponseCurve()
            response.gid = gid
            response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
            response.compound = row.cpd_name
            response.sample = row.ccl_name
            s = response.summary.add()
            s.type = phenotype_pb2.ResponseSummary.EC50
            s.value = row.apparent_ec50_umol
            s.unit = "uM"
            
            for m in data.loc[lambda x: x.master_cpd_id==row.master_cpd_id, : ].loc[lambda x: x.experiment_id==row.experiment_id].itertuples():
                dr = response.values.add()
                dr.dose = m.cpd_conc_umol
                dr.response = m.cpd_expt_avg_log2
            
            emit(response)
            gid_set.add(gid)


def convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, out, multi=None):


    # Read in Compound information into a pandas dataframe.
    compound_df = pandas.read_table(metadrugPath)
    # Read in Cell line information
    ccl_df = pandas.read_table(metacelllinePath)
    # Read in data curves for experiments
    datacurves_df = pandas.read_table(responsePath)
    # Read in meta experimental data
    metaexperiment_df = pandas.read_table(metaexperimentPath)
    
    ctdd_merged = pandas.merge(datacurves_df, metaexperiment_df, how='left', on=['experiment_id']) # merge experiment data
    ctdd_merged = pandas.merge(ctdd_merged, compound_df, how='left', on=['master_cpd_id']) # merge with compound data frame
    ctdd_merged = pandas.merge(ctdd_merged, ccl_df, how='left', on=['master_ccl_id']) # merge with cell line data frame
    
    ctdd_data = pandas.read_table(dataPath)
    #print ctdd_merged
    
    out_handles = {}
    def emit_json_single(message):
        if 'main' not in out_handles:
            out_handles['main'] = open(out, "w")
        msg = json.loads(json_format.MessageToJson(message))
        msg["#label"] = message.DESCRIPTOR.full_name
        out_handles['main'].write(json.dumps(msg))
        out_handles['main'].write("\n")
    def emit_json_multi(message):
        if message.DESCRIPTOR.full_name not in out_handles:
            out_handles[message.DESCRIPTOR.full_name] = open(multi + "." + message.DESCRIPTOR.full_name + ".json", "w")
        msg = json.loads(json_format.MessageToJson(message))
        out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
        out_handles[message.DESCRIPTOR.full_name].write("\n")
    if out is not None:
        emit = emit_json_single
    if multi is not None:
        emit = emit_json_multi
    
    ctdd_merged.to_csv("test.out", sep="\t")    
    
    process_drugs(emit, ctdd_merged)
    process_response(emit, ctdd_merged, ctdd_data)
    
########################################

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)


def convert_to_profobuf(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, out, multi):
    if responsePath and metadrugPath and metacelllinePath and metaexperimentPath and (out or multi) and format:
        convert_all_ctdd(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, dataPath, out, multi)
    else:
        print("Please include all arguments")

    #write_messages(state, outpath, format)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(options.response, options.metadrug, options.metacellline, options.metaexperiment, dataPath=options.data, out=options.out, multi=options.multi)

#! /usr/bin/python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

This program converts CTDD drug response information into
protobuf data based on the BMEG sample.proto schema.

Source: ftp://caftpd.nci.nih.gov/pub/dcc_ctd2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/ --> CTRPv2.0_2015_ctd2_ExpandedDataset.zip

The four files of interest (for this converter) are:
1) v20.data.curves_post_qc.txt
2) v20.meta.per_compound.txt
3) v20.meta.per_cell_line.txt
4) v20.meta.per_experiment.txt
'''

import sample_pb2 as schema
from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re
import pandas as pd

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
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

########################################

def find_biosample(state, source, barcode, sample_type):
    sample_name = 'biosample:' + source + '-' + barcode
    biosample = state['Biosample'].get(sample_name)
    if biosample is None:
        biosample = schema.Biosample()
        biosample.name = sample_name
        biosample.source = source
        biosample.barcode = barcode
        biosample.sampleType = sample_type
        state['Biosample'][sample_name] = biosample

    return biosample

def append_unique(l, i):
    if not i in l:
        l.append(i)

def process_row(state, source, row): #row is a namedtuple
    # format cell line name as in CCLE dataset
    ccl_name_CCLE = str(row.ccl_name) + "_" + str(row.ccle_primary_site).upper()
    sample_type = "tumor"
    # create Biosample message for cancer cell line
    biosample_source = "CCLE"
    tumor_sample = find_biosample(state, biosample_source, ccl_name_CCLE, sample_type)
    
    # Create OntologyTerm
    ontology_term_name = "ontologyTerm:" + "http://amigo.geneontology.org/amigo/term/GO:0042493"
    ontology_term = state['OntologyTerm'].get(ontology_term_name)
    if ontology_term is None:
        ontology_term = schema.OntologyTerm()
        ontology_term.name = ontology_term_name
        ontology_term.term = "response to drug"
        state['OntologyTerm'][ontology_term_name] = ontology_term
    # Create Phenotype based on the OntologyTerm
    phenotype_name = "phenotype:" + ontology_term.name #may want to refine this later
    phenotype = state['Phenotype'].get(phenotype_name)
    if phenotype is None:
        phenotype = schema.Phenotype()
        phenotype.name = phenotype_name
        state['Phenotype'][phenotype_name] = phenotype

    # create drug message for CTDD compound
    drug_name = "drug:" + row.cpd_name # in the future might want to find canonical drug name via external resources
    drug = state['Drug'].get(drug_name)
    if drug is None:
        drug = schema.Drug()
        drug.name = drug_name
        drug.synonyms.append("drug:" + row.cpd_name)
        drug.synonyms.append("drug:" + row.broad_cpd_id)
        drug.synonyms.append("drug:" + row.cpd_smiles)
        state['Drug'][drug_name] = drug

    # Create PhenotypeAssociation (which contains Context) based on the Drug
    phenotype_association_name = "phenotypeAssociation:" + tumor_sample.name + drug.name + phenotype.name
    phenotype_association = state['PhenotypeAssociation'].get(phenotype_association_name)
    if phenotype_association is None:
        phenotype_association = schema.PhenotypeAssociation()
        phenotype_association.name = phenotype_association_name
        phenotype_association.info['conc_pts_fit'] = str(row.conc_pts_fit)
        phenotype_association.info['fit_num_param'] = str(row.fit_num_param)
        phenotype_association.info['p1_conf_int_high'] = str(row.p1_conf_int_high)
        phenotype_association.info['p1_conf_int_low'] = str(row.p1_conf_int_low)
        phenotype_association.info['p2_conf_int_high'] = str(row.p2_conf_int_high)
        phenotype_association.info['p2_conf_int_low'] = str(row.p2_conf_int_low)
        phenotype_association.info['p4_conf_int_high'] = str(row.p4_conf_int_high)
        phenotype_association.info['p4_conf_int_low'] = str(row.p4_conf_int_low)
        phenotype_association.info['p1_center'] = str(row.p1_center)
        phenotype_association.info['p2_slope'] = str(row.p2_slope)
        phenotype_association.info['p3_total_decline'] = str(row.p3_total_decline)
        phenotype_association.info['p4_baseline'] = str(row.p4_baseline)
        phenotype_association.info['apparent_ec50_umol'] = str(row.apparent_ec50_umol)
        phenotype_association.info['pred_pv_high_conc'] = str(row.pred_pv_high_conc)
        phenotype_association.info['area_under_curve'] = str(row.area_under_curve)
        phenotype_association.info['run_id'] = str(row.run_id)
        phenotype_association.info['experiment_date'] = str(row.experiment_date)
        phenotype_association.info['culture_media'] = str(row.culture_media)
        phenotype_association.info['baseline_signal'] = str(row.baseline_signal)
        phenotype_association.info['cells_per_well'] = str(row.cells_per_well)
        phenotype_association.info['growth_mode'] = str(row.growth_mode)
        phenotype_association.info['snp_fp_status'] = str(row.snp_fp_status)
        phenotype_association.info['top_test_conc_umol'] = str(row.top_test_conc_umol)
        phenotype_association.info['cpd_status'] = str(row.cpd_status)
        phenotype_association.info['inclusion_rationale'] = str(row.inclusion_rationale)
        phenotype_association.info['gene_symbol_of_protein_target'] = str(row.gene_symbol_of_protein_target)
        phenotype_association.info['target_or_activity_of_compound'] = str(row.target_or_activity_of_compound)
        phenotype_association.info['source_name'] = str(row.source_name)
        phenotype_association.info['source_catalog_id'] = str(row.source_catalog_id)
        phenotype_association.info['ccl_availability'] = str(row.ccl_availability)
        phenotype_association.info['ccle_primary_hist'] = str(row.ccle_primary_hist)
        phenotype_association.info['ccle_hist_subtype_1'] = str(row.ccle_hist_subtype_1)
        state['PhenotypeAssociation'][phenotype_association_name] = phenotype_association

    # make edges
    append_unique(phenotype.isTypeEdges, ontology_term.name)
    append_unique(phenotype_association.hasGenotypeEdges, tumor_sample.name)
    append_unique(phenotype_association.hasPhenotypeEdges, phenotype.name)
    append_unique(phenotype_association.hasContextEdges, drug.name)

    return state

def convert_all_ctdd(state, source, responsePath, metadrugPath, metacelllinePath, metaexperimentPath):


    # Read in Compound information into a pandas dataframe.
    compound_df = pd.read_table(metadrugPath)
    # Read in Cell line information
    ccl_df = pd.read_table(metacelllinePath)
    # Read in data curves for experiments
    datacurves_df = pd.read_table(responsePath)
    # Read in meta experimental data
    metaexperiment_df = pd.read_table(metaexperimentPath)


    ctdd_merged = pd.merge(datacurves_df, metaexperiment_df, how='left', on=['experiment_id']) # merge experiment data
    ctdd_merged = pd.merge(ctdd_merged, compound_df, how='left', on=['master_cpd_id']) # merge with compound data frame
    ctdd_merged = pd.merge(ctdd_merged, ccl_df, how='left', on=['master_ccl_id']) # merge with cell line data frame

    # Iterate through each row of the merged dataframe and create protobuf messages when necessary
    for row in ctdd_merged.itertuples():
        process_row(state, source, row)

########################################

def splice_path(path, s):
    #print(path, s)
    path_split = path.split('.')
    suffix = path_split[-1]
    path_parts = path_split[:-1]
    path_parts.extend([s, suffix])
    return '.'.join(path_parts) #string.join(path_parts, '.')

def message_to_json(message):
    json = json_format.MessageToJson(message)
    return re.sub(r' +', ' ', json.replace('\n', ''))

def write_messages(state, outpath, format):
    if format == 'json':
        for message in state:
            outmessage = splice_path(outpath, message)
            messages = list(map(message_to_json, state[message].values())) #[message_to_json(value) for value in state[message].values()]
            if len(messages) > 0:
                out = '\n'.join(messages) #string.join(messages, '\n')
                outhandle = open(outmessage, 'w')
                outhandle.write(out)
                outhandle.close()
    else:
        for message in state:
            outmessage = splice_path(outpath, message)
            messages = list(map(lambda m: m.SerializeToString(), state[message].values())) #[value.SerializeToString() for value in state[message].values()]
            if len(messages) > 0:
                out = b'\n'.join(messages) #string.join(messages, '\n')
                outhandle = open(outmessage, 'wb')
                outhandle.write(out)
                outhandle.close()

def convert_to_profobuf(responsePath, metadrugPath, metacelllinePath, metaexperimentPath, outpath, format):
    state = {'Biosample': {},
             'Drug': {},
             'OntologyTerm': {},
             'Phenotype': {},
             'PhenotypeAssociation': {}}
    source = 'CTDD'

    if responsePath and metadrugPath and metacelllinePath and metaexperimentPath and outpath and format:
        convert_all_ctdd(state, source, responsePath, metadrugPath, metacelllinePath, metaexperimentPath)
    else:
        print("Please include all arguments")

    write_messages(state, outpath, format)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(options.response, options.metadrug, options.metacellline, options.metaexperiment, options.out, options.format)

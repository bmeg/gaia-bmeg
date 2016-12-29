#! /usr/bin/env python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

This program converts CCLE maf files and drug response information into
protobuf data based on the BMEG sample.proto schema.

Intended CCLE sources:
http://www.broadinstitute.org/ccle/data/browseData?conversationPropagation=begin#

Hybrid capture sequencing  (5.2MB) Mutation  (In process)  ->  CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz (preferred dataset)
Pharmacological profiling  (8.0MB) Drug data  (Published)  ->  CCLE_NP24.2009_Drug_data_2015.02.24.csv

CCLE Maf format contains some columns output by Oncotator: https://www.broadinstitute.org/oncotator/help/  (see Output Format section)


curl -o CCLE_NP24.2009_Drug_data_2015.02.24.csv "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_Drug_data_2015.02.24.csv?downloadff=true&fileId=20777"
curl -o CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_26/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz?downloadff=true&fileId=6873"

Clinical data download
curl -o CCLE_sample_info_file_2012-10-18.txt "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801"

Convert CCLE
curl -o CCLE_Expression_2012-09-29.res "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_21/CCLE_Expression_2012-09-29.res?downloadff=true&fileId=6760"

'''

from ga4gh import bio_metadata_pb2, variants_pb2, allele_annotations_pb2
from bmeg import phenotype_pb2, genome_pb2, matrix_pb2
from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re

########################################

def gid_biosample(name):
    return 'biosample:' + "CCLE:" + name

def gid_variant_set(name):
    return "variant_set:" + "CCLE:" + name

def gid_callset(name):
    return "callset:" + "CCLE:" + name

def gid_variant(chromosome, start, end, strand, ref, alt):
    return "variant:%s:%s:%s:%s:%s:%s" % ( chromosome, start, end, strand, ",".join(ref), ",".join(alt) )

def gid_compound(name):
    return "compound:" + name

def gid_response_curve(cellline, compound):
    return "responseCurve:%s:%s" % (cellline, compound)

def gid_expression(name):
    return "expression:%s" % (name)


########################################

def process_csv_line(state, source, lineAsList):
    # Fill out a PhenotypeAssociation connecting Biosample to a Drug response

    # Information indices
    CCLE_Cell_Line_Name = 0 #same as tumor_sample_barcode in process_maf_line()  #For Biosample
    Primary_Cell_Line_Name = 1
    Compound = 2 #for Drug
    # The rest of these fields go into the info field of the Context for this PhenotypeAssociation
    Target = 3
    Doses_uM = 4
    Activity_Data_median = 5
    Activity_SD = 6
    Num_Data = 7
    FitType = 8
    EC50_uM = 9
    IC50_uM = 10
    Amax = 11
    ActArea = 12

    # First create all the data structures to be referenced by the PhenotypeAssociation

    # Create Biosample
    tumor_sample = find_biosample(state, source, lineAsList[CCLE_Cell_Line_Name], 'cellline')

    # Create OntologyTerm
    """
    ontology_term_name = "ontologyTerm:" + "http://amigo.geneontology.org/amigo/term/GO:0042493"
    ontology_term = state['OntologyTerm'].get(ontology_term_name)
    if ontology_term is None:
        ontology_term = phenotype_pb2.OntologyTerm()
        ontology_term.gid = ontology_term_name
        ontology_term.type = 'OntologyTerm'
        ontology_term.term = "response to drug"
        state['OntologyTerm'][ontology_term_name] = ontology_term

    # Create Phenotype based on the OntologyTerm
    phenotype_name = "phenotype:" + ontology_term.gid #may want to refine this later
    phenotype = state['Phenotype'].get(phenotype_name)
    if phenotype is None:
        phenotype = phenotype_pb2.Phenotype()
        phenotype.gid = phenotype_name
        phenotype.type = 'Phenotype'
        state['Phenotype'][phenotype_name] = phenotype
    """
    

    """
    # Create PhenotypeAssociation (which contains Context) based on the Drug
    phenotype_association_name = "phenotypeAssociation:" + tumor_sample.gid + drug.gid + phenotype.gid
    phenotype_association = state['PhenotypeAssociation'].get(phenotype_association_name)
    if phenotype_association is None:
        phenotype_association = phenotype_pb2.PhenotypeAssociation()
        phenotype_association.gid = phenotype_association_name
        phenotype_association.type = 'PhenotypeAssociation'
        phenotype_association.infoProperties['Target'] = lineAsList[Target]
        phenotype_association.infoProperties['Doses_uM'] = lineAsList[Doses_uM]
        phenotype_association.infoProperties['Activity_Data_median'] = lineAsList[Activity_Data_median]
        phenotype_association.infoProperties['Activity_SD'] = lineAsList[Activity_SD]
        phenotype_association.infoProperties['Num_Data'] = lineAsList[Num_Data]
        phenotype_association.infoProperties['FitType'] = lineAsList[FitType]
        phenotype_association.infoProperties['EC50_uM'] = lineAsList[EC50_uM]
        phenotype_association.infoProperties['IC50_uM'] = lineAsList[IC50_uM]
        phenotype_association.infoProperties['Amax'] = lineAsList[Amax]
        phenotype_association.infoProperties['ActArea'] = lineAsList[ActArea]
        state['PhenotypeAssociation'][phenotype_association_name] = phenotype_association
    """

    return state

def convert_ccle_pharma_profiles(emit, csvpath):
    print('converting csv:' + csvpath)
    with open(csvpath) as pharma_file:
        
        drugs = {}
        
        reader = csv.DictReader(pharma_file, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
        for row in reader:
            drug_gid = gid_compound(row['Compound'])
            if drug_gid not in drugs:
                drug = phenotype_pb2.Compound()
                drug.gid = drug_gid
                drugs[drug_gid] = drug
            
            response = phenotype_pb2.ResponseCurve()
            response.gid = gid_response_curve(row['CCLE Cell Line Name'], row['Compound'])
            response.responseType = phenotype_pb2.ResponseCurve.ACTIVITY
            response.compound = drug_gid
            response.sample = gid_biosample(row['CCLE Cell Line Name'])
            for dose, activity in zip( row["Doses (uM)"].split(","), row["Activity Data (median)"].split(",") ):
                dr = response.values.add()
                dr.dose = float(dose)
                dr.response = float(activity)
            try:
                v = float(row["EC50 (uM)"])
                s = response.summary.add()
                s.type = phenotype_pb2.ResponseSummary.EC50
                s.value = v
                s.unit = "uM"
            except ValueError:
                pass
            try:
                v = float(row["IC50 (uM)"])
                s = response.summary.add()
                s.type = s.IC50
                s.value = v
                s.unit = "uM"
            except ValueError:
                pass
            emit(response)

        for d in drugs.values():
            emit(d)



def convert_expression(emit, expressionpath):

    with open(expressionpath) as handle:
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
            if len(line[0]):
                for a in zip(header, line)[2:]:
                    if len(a[0]):
                        values[a[0]].append(float(a[1]))
                gene_order.append(line[0])

        for k, v in values.items():
            if len(k):
                ge = matrix_pb2.GeneExpression()
                ge.gid = gid_expression(k)
                ge.bio_sample_id = gid_biosample(k)
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
                emit(ge)

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a
    

def convert_sample(emit, samplepath):
    
    with open(samplepath) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line in reader:
            sample = bio_metadata_pb2.BioSample()
            sample.id = gid_biosample(line["CCLE name"])
            sample.dataset_id = "CCLE"
            proto_list_append(sample.info['sampleType'], "cellline")
            proto_list_append(sample.info['histology'], line["Histology"])
            proto_list_append(sample.info['alias'], line["Cell line primary name"])
            if len(line['Source']):
                proto_list_append(sample.info['source'], line['Source'])
            if len(line['Notes']):
                proto_list_append(sample.info['notes'], line['Notes'])
            emit(sample)

########################################

def splice_path(path, s):
    #print(path, s)
    path_split = path.split('.')
    suffix = path_split[-1]
    path_parts = path_split[:-1]
    path_parts.extend([s, suffix])
    return '.'.join(path_parts) #string.join(path_parts, '.')

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)

def convert_to_profobuf(drugpath, samplepath, expressionpath, out, multi, format):
    
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
         
    if drugpath:
        convert_ccle_pharma_profiles(emit, drugpath)
    if samplepath:
        convert_sample(emit, samplepath)
    if expressionpath:
        convert_expression(emit, expressionpath)

    for handle in out_handles.values():
        handle.close()

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--drug', type=str, help='Path to the csv (drug response data) you want to import')
    parser.add_argument('--sample', type=str, help='Path to the csv CCLE sample data')
    parser.add_argument('--expression', type=str, help='Path to the CCLE expression data')
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(options.drug, options.sample, options.expression, options.out, options.multi, options.format)

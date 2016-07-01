#! /usr/bin/python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

This program converts CCLE maf files and drug response information into
protobuf data based on the BMEG sample.proto schema.

Intended CCLE sources:
http://www.broadinstitute.org/ccle/data/browseData?conversationPropagation=begin#

Hybrid capture sequencing  (6.5GB) Mutation  (In process)  ->  CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz (preferred dataset)
Pharmacological profiling  (8.0MB) Drug data  (Published)  ->  CCLE_NP24.2009_Drug_data_2015.02.24.csv

CCLE Maf format contains some columns output by Oncotator: https://www.broadinstitute.org/oncotator/help/  (see Output Format section)

'''

import sample_pb2 as schema
from google.protobuf import json_format
import json, sys, argparse, os
import csv #for drug data
import string
import re

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--maf', type=str, help='Path to the maf you want to import')
    parser.add_argument('--csv', type=str, help='Path to the csv (drug response data) you want to import')
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

def find_position(state, chromosome, start, end, strand):
    position_name = 'position:' + chromosome + start + end + strand
    position = state['Position'].get(position_name)
    if position is None:
        position = schema.Position()
        position.name = position_name
        position.chromosome = chromosome
        position.start = int(start)
        position.end = int(end)
        position.strand = strand
        state['Position'][position_name] = position

    return position

def find_feature(state, hugo_code, entrez_gene_id):
    feature_name = "feature:" + hugo_code
    feature = state['Feature'].get(feature_name)
    if feature is None:
        feature = schema.Feature()
        feature.name = feature_name
        feature.attributes['entrezGeneId'] = entrez_gene_id
        state['Feature'][feature_name] = feature

    return feature

def find_variant_call(state, source, position_name, reference_allele, normal_allele1, normal_allele2, tumor_allele1, tumor_allele2, variant_type, ncbi_build, mutation_status, sequencing_phase, sequence_source, bam_file, sequencer, genome_change, tumor_sample_name):
    variant_name = 'variantCall:' + source + tumor_sample_name + position_name + variant_type + mutation_status 
    variant_call = state['VariantCall'].get(variant_name)
    if variant_call is None:
        variant_call = schema.VariantCall()
        variant_call.name = variant_name
        variant_call.source = source
        
        variant_call.referenceAllele = reference_allele # Reference_Allele, e.g. 'T'
        variant_call.normalAllele1 = normal_allele1
        variant_call.normalAllele2 = normal_allele2
        variant_call.tumorAllele1 = tumor_allele1
        variant_call.tumorAllele2 = tumor_allele2
        variant_call.variantType = variant_type # e.g. 'SNP'
    
        variant_call.info['ncbiBuild'] = ncbi_build # e.g. '37'
        variant_call.info['mutationStatus'] = mutation_status # e.g. 'Somatic'
        variant_call.info['sequencingPhase'] = sequencing_phase # e.g. 'Phase_IV'
        variant_call.info['sequenceSource'] = sequence_source # e.g. 'Capture'
        variant_call.info['bamFile'] = bam_file # e.g. 'dbGAP'
        variant_call.info['sequencer'] = sequencer # e.g. 'dbGAP'
        variant_call.info['genomeChange'] = genome_change # e.g. 'dbGAP'
        state['VariantCall'][variant_name] = variant_call

    return variant_call

def find_variant_call_effect(state, source, effect_name, classification, line):

    dbsnp_rs = 13
    dbsnp_val_status = 14
    annotation_transcript = 33
    transcript_strand = 34
    cdna_change = 35 #not sure where this should belong
    codon_change = 36
    protein_change = 37
    other_transcripts = 38
    refseq_mrna_id = 39
    refseq_prot_id = 40
    swissprot_acc_id = 41
    swissprot_entry_id = 42
    description = 43
    uniprot_aapos = 44
    uniprot_region = 45
    uniprot_site = 46
    alternative_allele_reads_count = 47
    reference_allele_reads_count = 48
    x46vertebrates_aa_alignment_column = 49
    method = 50
 
    variant_call_effect = state['VariantCallEffect'].get(effect_name)
    if variant_call_effect is None:
        variant_call_effect = schema.VariantCallEffect()
        variant_call_effect.name = effect_name
        variant_call_effect.source = source
        variant_call_effect.variantClassification = classification # e.g. 'Missense_Mutation'
        try:
            variant_call_effect.info['dbsnp_rs'] = line[dbsnp_rs]
            variant_call_effect.info['dbsnp_val_status'] = line[dbsnp_val_status]
            variant_call_effect.info['annotation_transcript'] = line[annotation_transcript]
            variant_call_effect.info['transcript_strand'] = line[transcript_strand]
            variant_call_effect.info['cdna_change'] = line[cdna_change]
            variant_call_effect.info['codon_change'] = line[codon_change]
            variant_call_effect.info['protein_change'] = line[protein_change]
            variant_call_effect.info['other_transcripts'] = line[other_transcripts]
            variant_call_effect.info['refseq_mrna_id'] = line[refseq_mrna_id]
            variant_call_effect.info['refseq_prot_id'] = line[refseq_prot_id]
            variant_call_effect.info['swissprot_acc_id'] = line[swissprot_acc_id]
            variant_call_effect.info['swissprot_entry_id'] = line[swissprot_entry_id]
            variant_call_effect.info['description'] = line[description]
            variant_call_effect.info['uniprot_aapos'] = line[uniprot_aapos]
            variant_call_effect.info['uniprot_region'] = line[uniprot_region]
            variant_call_effect.info['uniprot_site'] = line[uniprot_site]
            variant_call_effect.info['alternative_allele_reads_count'] = line[alternative_allele_reads_count]
            variant_call_effect.info['reference_allele_reads_count'] = line[reference_allele_reads_count]
            variant_call_effect.info['x46vertebrates_aa_alignment_column'] = line[x46vertebrates_aa_alignment_column]
            variant_call_effect.info['method'] = line[method]
        except Exception as ex:
            None

        state['VariantCallEffect'][effect_name] = variant_call_effect


    #note to self: at this point, both variant_call_effect and state['VariantCallEffect'][effect_name] are bound to the same object (a protobuf message). mutating either object-bound-name will affect the final result.
    return variant_call_effect
    
def append_unique(l, i):
    if not i in l:
        l.append(i)

def process_maf_line(state, source, line_raw):
    # Information indices for VariantCall, Position, and Biosample
    center = 2
    ncbi_build = 3
    chromosome = 4
    start = 5
    end = 6
    strand = 7
    variant_type = 9
    reference_allele = 10
    tumor_allele1 = 11
    tumor_allele2 = 12
    tumor_sample_barcode = 15

    normal_sample_barcode = 16 #null for CCLE
    normal_allele1 = 17 #null for CCLE
    normal_allele2 = 18 #null for CCLE
    verification_status = 23 #null for CCLE
    validation_status = 24 #null for CCLE
    mutation_status = 25 #null for CCLE
    sequencing_phase = 26 #All Phase_I for CCLE
    sequence_source = 27 #Unspecified for CCLE
    bam_file = 30 #null for CCLE

    sequencer = 31
    genome_change = 32 #for varcal

    # Information indices for VariantCallEffect and Feature
    hugo_symbol = 0
    entrez_gene_id = 1
    variant_classification = 8

    # --------------------------------------------
    line = line_raw.rstrip().split('\t')

    # Example tumor_sample_barcode: 'CCLE-CCK81_LARGE_INTESTINE' 

    # create nodes
    tumor_sample = find_biosample(state, source, line[tumor_sample_barcode], 'tumor')

    position = find_position(state, line[chromosome], line[start], line[end], line[strand])

    feature = find_feature(state, line[hugo_symbol], line[entrez_gene_id])

    variant_call = find_variant_call(state, source, position.name, line[reference_allele], line[normal_allele1], line[normal_allele2], line[tumor_allele1], line[tumor_allele2], line[variant_type], line[ncbi_build], line[mutation_status], line[sequencing_phase], line[sequence_source], line[bam_file], line[sequencer], line[genome_change], tumor_sample.name)

    effect_name = "variantCallEffect:" + variant_call.name #For now, each variant call has one variant call effect. In the future we might want to make the effect name more unique.
    variant_call_effect = find_variant_call_effect(state, source, effect_name, line[variant_classification], line)

    # make edges
    append_unique(variant_call.atPositionEdges, position.name)
    append_unique(variant_call.tumorSampleEdges, tumor_sample.name)

    append_unique(variant_call_effect.inFeatureEdges, feature.name)
    append_unique(variant_call_effect.effectOfEdges, variant_call.name)

    return state

def convert_maf(state, mafpath, source):
    print('converting maf: ' + mafpath)

    inhandle = open(mafpath)
    for line in inhandle:
        if not line.startswith('Hugo_Symbol') and not line.startswith('#'):
            state = process_maf_line(state, source, line)
    inhandle.close()

    return state

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
    tumor_sample = find_biosample(state, source, lineAsList[CCLE_Cell_Line_Name], 'tumor')

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

    # Create Drug
    drug_name = "drug:" + lineAsList[Compound]
    drug = state['Drug'].get(drug_name)
    if drug is None:
        drug = schema.Drug()
        drug.name = drug_name
        drug.synonyms.append(drug_name)
        state['Drug'][drug_name] = drug
    # Create PhenotypeAssociation (which contains Context) based on the Drug
    phenotype_association_name = "phenotypeAssociation:" + tumor_sample.name + drug.name + phenotype.name
    phenotype_association = state['PhenotypeAssociation'].get(phenotype_association_name)
    if phenotype_association is None:
        phenotype_association = schema.PhenotypeAssociation()
        phenotype_association.name = phenotype_association_name
        phenotype_association.info['Target'] = lineAsList[Target]
        phenotype_association.info['Doses_uM'] = lineAsList[Doses_uM]
        phenotype_association.info['Activity_Data_median'] = lineAsList[Activity_Data_median]
        phenotype_association.info['Activity_SD'] = lineAsList[Activity_SD]
        phenotype_association.info['Num_Data'] = lineAsList[Num_Data]
        phenotype_association.info['FitType'] = lineAsList[FitType]
        phenotype_association.info['EC50_uM'] = lineAsList[EC50_uM]
        phenotype_association.info['IC50_uM'] = lineAsList[IC50_uM]
        phenotype_association.info['Amax'] = lineAsList[Amax]
        phenotype_association.info['ActArea'] = lineAsList[ActArea]
        state['PhenotypeAssociation'][phenotype_association_name] = phenotype_association

    # make edges
    append_unique(phenotype.isTypeEdges, ontology_term.name)
    append_unique(phenotype_association.hasGenotypeEdges, tumor_sample.name)
    append_unique(phenotype_association.hasPhenotypeEdges, phenotype.name)
    append_unique(phenotype_association.hasContextEdges, drug.name)

    return state

def convert_ccle_pharma_profiles(state, csvpath, source):
    print('converting csv:' + csvpath)
    with open(csvpath) as pharma_file:
        reader = csv.reader(pharma_file, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
        for row in reader:
            if not row[0] == 'CCLE Cell Line Name':
                state = process_csv_line(state, source, row)

    return state

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

def convert_maf_and_csv_to_profobuf(mafpath, csvpath, outpath, format):
    state = {'Biosample': {},
             'Position': {},
             'Feature': {},
             'VariantCall': {},
             'VariantCallEffect': {},
             'Drug': {},
             'OntologyTerm': {},
             'Phenotype': {},
             'PhenotypeAssociation': {}}
    source = 'CCLE'

    if mafpath:
        convert_maf(state, mafpath, source)
    if csvpath:
        convert_ccle_pharma_profiles(state, csvpath, source)

    write_messages(state, outpath, format)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_maf_and_csv_to_profobuf(options.maf, options.csv, options.out, options.format)

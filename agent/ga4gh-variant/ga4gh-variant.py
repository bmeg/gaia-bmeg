#! /usr/bin/env python
'''
Authors: Malisa Smith smimal@ohsu.edu, Ryan Spangler spanglry@ohsu.edu

'''

from ga4gh import bio_metadata_pb2, variants_pb2, allele_annotations_pb2
from google.protobuf import json_format
import json, sys, argparse, os
import string
import logging
import re

########################################

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a


class Converter(object):
    def __init__(self, bioPrefix, variantPrefix, variantSetPrefix, callsetPrefix, variantAnnPrefix, transEffectPrefix, hugoPrefix, typeField):
        self.bioPrefix = bioPrefix
        self.variantPrefix = variantPrefix
        self.variantSetPrefix = variantSetPrefix
        self.callsetPrefix = callsetPrefix
        self.typeField = typeField
        self.variantAnnPrefix = variantAnnPrefix
        self.transEffectPrefix = transEffectPrefix
        self.hugoPrefix = hugoPrefix
        
    def gid_biosample(self, name):
        return self.bioPrefix + name

    def gid_variant_set(self, name):
        return self.variantSetPrefix + name

    def gid_callset(self, name):
        return self.callsetPrefix + name

    def gid_variant(self, chromosome, start, end, strand, ref, alt):
        return "%s%s:%s:%s:%s:%s:%s" % ( self.variantPrefix, chromosome, start, end, strand, ",".join(ref), ",".join(alt) )
        
    def gid_variant_annotation(self, chromosome, start, end, strand, ref, alt):
        return "%s%s:%s:%s:%s:%s:%s" % ( self.variantAnnPrefix, chromosome, start, end, strand, ",".join(ref), ",".join(alt) )

    def gid_trans_effect(self, chromosome, start, end, strand, ref, alt):
        return "%s%s:%s:%s:%s:%s:%s" % ( self.transEffectPrefix, chromosome, start, end, strand, ",".join(ref), ",".join(alt) )

    def gid_gene(self, name):
        return "%s%s" % (self.hugoPrefix, name)

class MafConverter(Converter):
    #def __init__(self, **kwargs):
    #    super(MafConverter, self).__init__(**kwargs)

    def convert(self, emit, mafpath):
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

        normal_sample_barcode = 16 
        normal_allele1 = 17 
        normal_allele2 = 18 
        verification_status = 23
        validation_status = 24
        mutation_status = 25
        sequencing_phase = 26
        sequence_source = 27
        bam_file = 30

        sequencer = 31
        genome_change = 32

        # Information indices for VariantCallEffect and Gene
        hugo_symbol = 0
        entrez_gene_id = 1
        variant_classification = 8
        
        logging.info('converting maf: ' + mafpath)
        
        samples = set()
        variants = {}
        variant_anns = {}

        inhandle = open(mafpath)
        for line_raw in inhandle:
            if not line_raw.startswith('Hugo_Symbol') and not line_raw.startswith('#'):
                line = line_raw.rstrip().split('\t')
                vid = self.gid_variant( line[chromosome], long(line[start]), long(line[end]), line[strand], 
                    line[reference_allele], set( [line[tumor_allele1], line[tumor_allele2] ] ) )
                if vid not in variants:
                    variant = variants_pb2.Variant()
                    variant.id = vid
                    variant.start = long(line[start])
                    variant.end = long(line[end])
                    variant.reference_name = line[chromosome]
                    variant.reference_bases = line[reference_allele]
                    for a in set( [line[tumor_allele1], line[tumor_allele2] ] ):
                        variant.alternate_bases.append(a)
                    proto_list_append(variant.info["hugo"], line[hugo_symbol])
                    proto_list_append(variant.info["center"], line[center])
                    proto_list_append(variant.info["variant_type"], line[variant_type])
                    variants[vid] = variant
                    
                    var_ann = allele_annotations_pb2.VariantAnnotation()
                    var_ann.id = self.gid_variant_annotation(line[chromosome], long(line[start]), long(line[end]), line[strand], 
                        line[reference_allele], set( [line[tumor_allele1], line[tumor_allele2] ] ) )
                    var_ann.variant_id = vid
                    trans_effect = var_ann.transcript_effects.add()
                    trans_effect.id = self.gid_trans_effect( line[chromosome], long(line[start]), long(line[end]), line[strand], 
                        line[reference_allele], set( [line[tumor_allele1], line[tumor_allele2] ] ) )
                    eff = trans_effect.effects.add()
                    eff.term = line[variant_type]
                    trans_effect.feature_id = self.gid_gene(line[hugo_symbol])
                    
                    variant_anns[vid] = var_ann
                    
                call = variants[vid].calls.add()
                call.call_set_id = self.gid_callset(line[tumor_sample_barcode])                
                samples.add(line[tumor_sample_barcode])

        for variant in variants.values():
            emit(variant)

        for variant in variant_anns.values():
            emit(variant)
        
        for i in samples:
            callset = variants_pb2.CallSet()
            callset.id = self.gid_callset(i)
            callset.biosample_id = self.gid_biosample(i)
            emit(callset)

        inhandle.close()


def convert_to_profobuf(maf, vcf, out, multi, format, bioPrefix, variantPrefix, variantSetPrefix, callsetPrefix, variantAnnPrefix, transEffectPrefix, hugoPrefix, typeField):
    out_handles = {}
    if maf:
        m = MafConverter(
            bioPrefix=bioPrefix,
            variantPrefix=variantPrefix,
            variantSetPrefix=variantSetPrefix,
            callsetPrefix=callsetPrefix,
            variantAnnPrefix=variantAnnPrefix,
            transEffectPrefix=transEffectPrefix,
            hugoPrefix=hugoPrefix,
            typeField=typeField
        )
        def emit_json_single(message):
            if 'main' not in out_handles:
                out_handles['main'] = open(out, "w")
            msg = json.loads(json_format.MessageToJson(message))
            msg[typeField] = message.DESCRIPTOR.full_name
            out_handles['main'].write(json.dumps(msg))
            out_handles['main'].write("\n")
        def emit_json_multi(message):
            if message.DESCRIPTOR.full_name not in out_handles:
                out_handles[message.DESCRIPTOR.full_name] = open(multi + "." + message.DESCRIPTOR.full_name + ".json", "w")
            msg = json.loads(json_format.MessageToJson(message))
            out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
            out_handles[message.DESCRIPTOR.full_name].write("\n")
        if out is not None:
            m.convert(emit_json_single, maf)
        if multi is not None:
            m.convert(emit_json_multi, maf)
    for handle in out_handles.values():
        handle.close()

def parse_args(args):
    # We don't need the first argument, which is the program name
    args = args[1:]

    # Construct the parser
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Now add all the options to it
    parser.add_argument('--maf', type=str, help='Path to the maf you want to import')
    parser.add_argument('--vcf', type=str, help='Path to the vcf you want to import')
    parser.add_argument('--out', type=str, help='Path to output file (.json or .pbf_ext)')
    parser.add_argument('--multi', type=str, help='Path to output file prefix (.json or .pbf_ext)')
    parser.add_argument('--bioPrefix', default="sample:")
    parser.add_argument('--variantPrefix', default="variant:")
    parser.add_argument('--variantSetPrefix', default="variantSet:")
    parser.add_argument('--callsetPrefix', default="callset:")
    parser.add_argument('--typeField', default="#label")
    parser.add_argument('--variantAnnPrefix', default="variantAnnotation:")
    parser.add_argument('--transEffectPrefix', default="transEffect:")
    parser.add_argument('--hugoPrefix', default="gene:")
    
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(**vars(options))

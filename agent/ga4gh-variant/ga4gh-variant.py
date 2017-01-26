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
import gzip

########################################

def proto_list_append(message, a):
    v = message.values.add()
    v.string_value = a


class Converter(object):
    def __init__(self, bioPrefix, variantPrefix, variantSetPrefix, callSetPrefix, variantAnnotationPrefix, transcriptEffectPrefix, hugoPrefix, typeField, centerCol):
        self.bioPrefix = bioPrefix
        self.variantPrefix = variantPrefix
        self.variantSetPrefix = variantSetPrefix
        self.callSetPrefix = callSetPrefix
        self.typeField = typeField
        self.variantAnnotationPrefix = variantAnnotationPrefix
        self.transcriptEffectPrefix = transcriptEffectPrefix
        self.hugoPrefix = hugoPrefix
        self.centerCol = centerCol
        
    def gid_biosample(self, name):
        return '%s:%s' % (self.bioPrefix, name)

    def gid_variant_set(self, name):
        return '%s:%s' % (self.variantSetPrefix, name)

    def gid_callSet(self, name):
        return '%s:%s' % (self.callSetPrefix, name)

    def gid_variant(self, chromosome, start, end, reference, alternate):
        return '%s:%s:%s:%s:%s:%s' % (self.variantPrefix, chromosome, start, end, reference, alternate)
        
    def gid_variant_annotation(self, variant_id, annotation_set_id):
        return '%s:%s:%s' % (self.variantAnnotationPrefix, variant_id, annotation_set_id)

    def gid_transcript_effect(self, feature_id, annotation_id, alternate_bases):
        return '%s:%s:%s:%s' % (self.transcriptEffectPrefix, feature_id, annotation_id, alternate_bases)

    def gid_gene(self, name):
        return '%s:%s' % (self.hugoPrefix, name)

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
        variants = set()

        header = None
        if mafpath.endswith(".gz"):
            inhandle = gzip.GzipFile(mafpath)
        else:
            inhandle = open(mafpath)
        for line_raw in inhandle:
            if line_raw.startswith('Hugo_Symbol'):
                header = line_raw.rstrip().split("\t")
            elif line_raw.startswith('#'):
                pass
            else:
                line = line_raw.rstrip().split('\t')
                row = dict(zip(header, line))
                alternate_bases = set([line[tumor_allele1], line[tumor_allele2]])
                variant_id = self.gid_variant(
                    line[chromosome],
                    long(line[start]),
                    long(line[end]),
                    line[reference_allele],
                    ','.join(alternate_bases))


                variant = variants_pb2.Variant()
                variant.id = variant_id
                variant.start = long(line[start])
                variant.end = long(line[end])
                variant.reference_name = line[chromosome]
                variant.reference_bases = line[reference_allele]
                for a in set( [line[tumor_allele1], line[tumor_allele2] ] ):
                    variant.alternate_bases.append(a)
                #proto_list_append(variant.info["hugo"], line[hugo_symbol])
                for c in row[self.centerCol].split("|"):
                    proto_list_append(variant.info["center"], c)
                #proto_list_append(variant.info["variant_type"], line[variant_type])
                call = variant.calls.add()
                call.call_set_id = self.gid_callSet(line[tumor_sample_barcode])
                samples.add(line[tumor_sample_barcode])
                emit(variant)

                annotation_id = self.gid_variant_annotation(variant_id, '')
                annotation = allele_annotations_pb2.VariantAnnotation()
                annotation.id = annotation_id
                annotation.variant_id = variant_id
                feature_id = self.gid_gene(line[hugo_symbol])
                effect = line[variant_type]
                transcript_effect = annotation.transcript_effects.add()
                transcript_effect.alternate_bases = ','.join(alternate_bases)
                transcript_effect.feature_id = feature_id
                transcript_effect.id = self.gid_transcript_effect(
                    feature_id,
                    annotation_id,
                    ','.join(alternate_bases))

                ontology = transcript_effect.effects.add()
                ontology.term = effect
                #annotations.append(annotation)
                emit(annotation)
        
        for i in samples:
            callset = variants_pb2.CallSet()
            callset.id = self.gid_callset(i)
            callset.bio_sample_id = self.gid_biosample(i)
            emit(callset)

        inhandle.close()


def convert_to_profobuf(maf, vcf, out, multi, format, bioPrefix, variantPrefix, variantSetPrefix, callSetPrefix, variantAnnotationPrefix, transcriptEffectPrefix, hugoPrefix, typeField, centerCol):
    out_handles = {}
    if maf:
        m = MafConverter(
            bioPrefix=bioPrefix,
            variantPrefix=variantPrefix,
            variantSetPrefix=variantSetPrefix,
            callSetPrefix=callSetPrefix,
            variantAnnotationPrefix=variantAnnotationPrefix,
            transcriptEffectPrefix=transcriptEffectPrefix,
            hugoPrefix=hugoPrefix,
            typeField=typeField,
            centerCol=centerCol
        )
        def emit_json_single(message):
            if 'main' not in out_handles:
                out_handles['main'] = open(out, 'w')
            msg = json.loads(json_format.MessageToJson(message))
            msg[typeField] = message.DESCRIPTOR.full_name
            out_handles['main'].write(json.dumps(msg))
            out_handles['main'].write("\n")
        def emit_json_multi(message):
            if message.DESCRIPTOR.full_name not in out_handles:
                out_handles[message.DESCRIPTOR.full_name] = open(multi + '.' + message.DESCRIPTOR.full_name + '.json', 'w')
            msg = json.loads(json_format.MessageToJson(message))
            out_handles[message.DESCRIPTOR.full_name].write(json.dumps(msg))
            out_handles[message.DESCRIPTOR.full_name].write('\n')
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
    parser.add_argument('--bioPrefix', default='biosample')
    parser.add_argument('--variantPrefix', default='variant')
    parser.add_argument('--variantSetPrefix', default='variantSet')
    parser.add_argument('--callSetPrefix', default='callSet')
    parser.add_argument('--typeField', default='#label')
    parser.add_argument('--variantAnnotationPrefix', default='variantAnnotation')
    parser.add_argument('--transcriptEffectPrefix', default='transcriptEffect')
    parser.add_argument('--hugoPrefix', default='gene')
    parser.add_argument('--center', dest="centerCol", default='center', help='caller field')
    
    parser.add_argument('--format', type=str, default='json', help='Format of output: json or pbf (binary)')
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_to_profobuf(**vars(options))

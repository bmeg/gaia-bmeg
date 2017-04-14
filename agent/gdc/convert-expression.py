#!/usr/bin/env python

import os
import re
import sys
import csv
import gzip
import json
import string
import argparse
from sets import Set
from pprint import pprint
from google.protobuf import json_format

import ga4gh.bio_metadata_pb2
import bmeg.matrix_pb2
import gdc_scan as scan
from google.protobuf import json_format

# example usage:
# python -m convert.gdc.convert-expression --ensembl ~/Data/hugo/ensembl.json --path ~/Data/gdc/prad/source/ --out ~/Data/gdc/prad/schema/prad-expression.json --tree ~/Data/gdc/gene-expression-file-samples.json

class ExpressionGenerator:
    def __init__(self, sample_gid):
        self.sample_gid = sample_gid

    def schema(self):
        return bmeg.matrix_pb2.GeneExpression()

    # This may need to be more specific (including aliquot?)
    def gid(self, data):
        return 'geneExpression:' + data['sample']['submitter_id']

    def update(self, expression, data):
        expression.barcode = data['sample']['submitter_id']
        for symbol in data['expression']:
            expression.expressions[symbol] = data['expression'][symbol]

        record.append_unique(expression.expressionForEdges, self.sample_gid(data))
        return expression

def file_tree():
    args = lambda: None
    args.type = 'Gene Expression Quantification'
    args.legacy = False
    tree = scan.case_files(args)
    return tree

def ensembl_hugo(path):
    with open(path) as ensembl:
        return json.loads(ensembl.read())

def translate_values(ensembl, raw):
    expression = {}

    symbols = []
    for line in raw.split('\n'):
        parts = line.split('\t')
        if len(parts) > 1:
            full_symbol, value = parts
            symbol = full_symbol.split('.')[0]
            symbols.append(symbol)
            if symbol in ensembl:
                expression[ensembl[symbol]] = float(value)
            else:
                pass # print('symbol not found: ' + symbol)

    return expression

def process_expression(tree, ensembl, file, raw, emit):
    print('processing ' + file)

    if file in tree:
        sample = tree[file]['samples'][0]
        expression = translate_values(ensembl, raw)
        
        out = bmeg.matrix_pb2.GeneExpression()
        for k,v in expression.items():
            if v != 0.0:
                out.expressions[k] = v
        out.biosample_id = 'biosample:' + tree[file]['project_id'] + ':' + sample['submitter_id']
        out.type = sample['sample_type']
        emit(out)
        

def convert_expression(path, ensembl, tree, emit):
    """
    def sample_gid(data):
        return 'biosample:' + data['sample']['submitter_id']

    state = {
        'GeneExpression': {},
        'types': ['GeneExpression'],
        'generators': {
            'expression': ExpressionGenerator(sample_gid)
        }
    }
    """
    
    files = os.listdir(path)
    failed = []

    print('iterating through files')
    for file in files:
        if re.search(r'\.FPKM\.', file):
            opener = open
            if file[-3:] == '.gz':
                opener = gzip.open

            try:
                with opener(os.path.join(path, file), 'rb') as f:
                    process_expression(tree, ensembl, file, f.read(), emit)
            except:
               failed.append(file)

    print('failed!\n' + str(failed))

    return state

def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)

def convert(options):
    print('fetching tree')
    if options.tree == 'gdc':
        tree = file_tree()
    else:
        with open(options.tree) as tree_file:
            tree = json.loads(tree_file.read())

    print('mapping ensembl to hugo')
    ensembl = ensembl_hugo(options.ensembl)
    handle = open(options.out, "w")
    def emit(record):
        handle.write(message_to_json(record) + "\n")
    convert_expression(options.path, ensembl, tree, emit)
    handle.close()

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ensembl', type=str, help='path to ensembl to hugo mapping')
    parser.add_argument('--path', type=str, help='path to expression files')
    parser.add_argument('--out', type=str, help='path to output file')
    parser.add_argument('--tree', type=str, default='gdc', help='path to case tree')
    
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert(options)

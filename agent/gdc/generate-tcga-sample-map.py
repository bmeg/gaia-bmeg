import os
import json
import argparse

def process_input(input, output):
    mapping = {}
    with open(input) as inhandle:
        for line in inhandle:
            sample = json.loads(line)
            mapping[sample['name']] = sample['datasetId']
    out = json.dumps(mapping)
    with open(output, 'w') as outhandle:
        outhandle.write(out)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='path to file containing tcga biosample messages')
    parser.add_argument('--output', type=str, help='path to output file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    process_input(args.input, args.output)

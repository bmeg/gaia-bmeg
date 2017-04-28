import re
import os
import sys
import json
import argparse

import bmeg.matrix_pb2 as schema
from google.protobuf import json_format

# import record as record

# fetch data ------------
# python gdc-scan.py cases list

def append_unique(l, i):
    if not i in l:
        l.append(i)

def initial_state(generators):
    state = {}
    state['types'] = generators.keys()
    state['generators'] = generators
    for key in generators:
        state[key] = {}

    return state

def message_to_json(message):
    json = json_format.MessageToJson(message)
    return re.sub(r' +', ' ', json.replace('\n', ''))

def process_input(state, path, process):
    if os.path.isdir(path):
        for item in os.listdir(path):
            with open(item) as input:
                state = process(state, input)
    else:
        with open(path) as input:
            state = process(state, input)

    return state

def output_state(state, path):
    json = []
    for type in state['types']:
        for key in state[type]:
            message = message_to_json(state[type][key])
            json.append(message)
    out = '\n'.join(json)
    outhandle = open(path, 'w')
    outhandle.write(out)
    outhandle.close()

class RecordGenerator(object):
    def __init__(self, name):
        self.name = name

    def schema(self):
        raise Exception('schema() not implemented')

    def gid(self, data):
        raise Exception('gid() not implemented')

    def update(self, record, data):
        raise Exception('update() not implemented')

    def create(self, data):
        record = self.schema()
        gid = self.gid(data)
        record.gid = gid
        record.type = self.name
        self.update(record, data)
        return record

    def find(self, state, data):
        gid = self.gid(data)
        record = state[self.name].get(gid)
        if record is None:
            record = self.create(data)
            state[self.name][gid] = record
        else:
            self.update(record, data)
        return record

def sample_gid(project, sample):
    return 'biosample:' + project + ":" + sample

class CohortGenerator(RecordGenerator):
    def __init__(self):
        super(CohortGenerator, self).__init__('Cohort')

    def schema(self):
        return schema.Cohort()

    def gid(self, data):
        return 'cohort:' + data['project']['project_id']

    def update(self, cohort, data):
        if not cohort.name:
            cohort.name = data['project']['project_id']
        if not cohort.location:
            cohort.location = data['project']['primary_site']
        if not cohort.description:
            cohort.description = data['project']['name']

        if 'submitter_sample_ids' in data:
            for sample in data['submitter_sample_ids']:
                append_unique(cohort.hasSample, sample_gid(cohort.name, sample))

        return cohort

def read_input(path):
    with open(path) as input:
        return input.read().split('\n')

def convert_cohort(input, output):
    state = {
        'Cohort': {},
        'types': ['Cohort'],
        'generators': {
            'cohort': CohortGenerator()
        }
    }

    lines = read_input(input)
    for line in lines:
        if not line == '':
            try:
                sample = json.loads(line)
            except:
                print('sample failed! ' + line)

            # print(sample)
            state['generators']['cohort'].find(state, sample)

    output_state(state, output)

def parse_args(args):
    args = args[1:]
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=str, help='path to input file')
    parser.add_argument('--output', type=str, help='path to output file')
    
    return parser.parse_args(args)

if __name__ == '__main__':
    options = parse_args(sys.argv)
    convert_cohort(options.input, options.output)

#!/usr/bin/env python

# download:
# ../gdc_scan.py files download --project TCGA-LUAD --type "Clinical Supplement"
# example usage:
# python -m convert.gdc.convert-clinical Individual --input ~/Data/gdc/clinical/ --output ~/Data/gdc/individual.json

import os
import re
import json
from copy import deepcopy
import argparse
from xml.dom.minidom import parseString
from google.protobuf import json_format

from ga4gh import bio_metadata_pb2

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

def dom_scan(node, query):
    stack = query.split('/')
    if node.localName == stack[0]:
        return dom_scan_iter(node, stack[1:], [stack[0]])

def dom_scan_iter(node, stack, prefix):
    if len(stack):
        for child in node.childNodes:
                if child.nodeType == child.ELEMENT_NODE:
                    if child.localName == stack[0]:
                        for out in dom_scan_iter(child, stack[1:], prefix + [stack[0]]):
                            yield out
                    elif '*' == stack[0]:
                        for out in dom_scan_iter(child, stack[1:], prefix + [child.localName]):
                            yield out
    else:
        if node.nodeType == node.ELEMENT_NODE:
            yield node, prefix, dict(node.attributes.items()), getText( node.childNodes )
        elif node.nodeType == node.TEXT_NODE:
            yield node, prefix, None, getText( node.childNodes )



def record_initial_state(generators):
    state = {}
    state['types'] = generators.keys()
    state['generators'] = generators
    for key in generators:
        state[key] = {}

    return state


def process_input(state, path, process):
    if os.path.isdir(path):
        for item in os.listdir(path):
            with open(os.path.join(path, item)) as input:
                state = process(state, input)
    else:
        with open(path) as input:
            state = process(state, input)

    return state


def message_to_json(message):
    msg = json.loads(json_format.MessageToJson(message))
    msg['#label'] = message.DESCRIPTOR.name
    return json.dumps(msg)

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
        record.id = gid
        #record.type = self.name
        self.update(record, data)
        return record

    def find(self, data):
        gid = self.gid(data)
        record = state[self.name].get(gid)
        if record is None:
            record = self.create(data)
            state[self.name][gid] = record
        else:
            self.update(record, data)
        return record

class IndividualGenerator(RecordGenerator):
    def __init__(self):
        super(IndividualGenerator, self).__init__('Individual')
        
    def schema(self):
        return bio_metadata_pb2.Individual()

    def gid(self, data):
        return 'individual:' + data['project_code'] + '-' + data['disease_code'] + ':' + data['bcr_patient_barcode']

    def update(self, individual, data):
        individual.name = data['bcr_patient_barcode']
        individual.dataset_id = data['project_code'] + '-' + data['disease_code']
        if 'submitted_tumor_site' in data:
            individual.info['tumorSite'].append(data['submitted_tumor_site'])

        for key in data:
            if len(data[key]):
                individual.info[key].append(data[key])
        return individual

class BiosampleGenerator(RecordGenerator):
    def __init__(self, individual_gid):
        super(BiosampleGenerator, self).__init__('Biosample')
        self.individual_gid = individual_gid

    def schema(self):
        return bio_metadata_pb2.Biosample()

    def gid(self, data):
        return 'biosample:' + data['project_code'] + '-' + data['disease_code'] + ':' + data['bcr_sample_barcode']

    def update(self, sample, data):
        sample.name = data['bcr_sample_barcode']
        sample.dataset_id = data['project_code'] + '-' + data['disease_code']
        if 'submitted_tumor_site' in data:
            site = data['submitted_tumor_site']
            sample.disease.term = site

        # sample.sampleType = 'tumor' if re.search('tumor', data['sample_type'], re.IGNORECASE) else 'normal'
        for key in data:
            if len(data[key]):
                sample.info[key].append(data[key])

        individual_id = {
            'project_code': data['project_code'],
            'disease_code': data['disease_code'],
            'bcr_patient_barcode': data['bcr_sample_barcode'][:12]
        }

        sample.individual_id = self.individual_gid(individual_id)
        # record.append_unique(sample.sampleOf, self.individual_gid({
        #     'bcr_patient_barcode': data['bcr_sample_barcode'][:12]}))

        return sample

def extract_attribute(data, stack, attr, text):
    if 'xsd_ver' in attr:
        p_name = attr.get('preferred_name', stack[-1])
        if len(p_name) == 0:
            p_name = stack[-1]
        data[p_name] = text

class ClinicalParser:
    def __init__(self):
        pass
    
    def parseXMLFile(self, state, dom, subtype):
        root_node = dom.childNodes[0]
        admin = {}
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/admin/*'):
            admin[stack[-1]] = text
        
        patient_barcode = None
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/bcr_patient_barcode'):
            patient_barcode = text
        
        patient_data = admin
        for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/*'):
            extract_attribute(patient_data, stack, attr, text)

        if subtype == 'Individual':
            for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/*'):
                extract_attribute(patient_data, stack, attr, text)
            for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/*/*'):
                extract_attribute(patient_data, stack, attr, text)
            for node, stack, attr, text in dom_scan(root_node, 'tcga_bcr/patient/stage_event/tnm_categories/*/*'):
                extract_attribute(patient_data, stack, attr, text)
            self.emit( state, patient_barcode, patient_data, subtype )
        
        if subtype == 'Biosample':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample'):
                sample_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'sample/bcr_sample_barcode'):
                    sample_barcode = c_text
                sample_data = deepcopy(patient_data)
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'sample/*'):
                    if 'xsd_ver' in c_attr:
                        sample_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, sample_barcode, sample_data, subtype )

        if subtype == 'portion':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion'):
                portion_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'portion/bcr_portion_barcode'):
                    portion_barcode = c_text
                portion_data = {}    
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'portion/*'):
                    if 'xsd_ver' in c_attr:
                        portion_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, portion_barcode, portion_data, 'Portion' )
        
        if subtype == 'analyte':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte'):
                analyte_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'analyte/bcr_analyte_barcode'):
                    analyte_barcode = c_text
                analyte_data = {}    
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'analyte/*'):
                    if 'xsd_ver' in c_attr:
                        analyte_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, analyte_barcode, analyte_data, 'Analyte' )

        if subtype == 'aliquot':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/samples/sample/portions/portion/analytes/analyte/aliquots/aliquot'):
                aliquot_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'aliquot/bcr_aliquot_barcode'):
                    aliquot_barcode = c_text
                aliquot_data = {}    
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'aliquot/*'):
                    if 'xsd_ver' in c_attr:
                        aliquot_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, aliquot_barcode, aliquot_data, 'Aliquot' )
        
        if subtype == 'drug':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/drugs/drug'):
                drug_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'drug/bcr_drug_barcode'):
                    drug_barcode = c_text
                drug_data = {}    
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'drug/*'):
                    if 'xsd_ver' in c_attr:
                        drug_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, drug_barcode, drug_data, 'Drug' )

        if subtype == 'radiation':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/radiations/radiation'):
                radiation_barcode = None
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'radiation/bcr_radiation_barcode'):
                    radiation_barcode = c_text
                radiation_data = {}    
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'radiation/*'):
                    if 'xsd_ver' in c_attr:
                        radiation_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, radiation_barcode, radiation_data, 'Radiation' )

        if subtype == 'followup':
            for s_node, s_stack, s_attr, s_text in dom_scan(root_node, 'tcga_bcr/patient/follow_ups/follow_up'):
                follow_up_barcode = None
                sequence = s_attr['sequence']
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'follow_up/bcr_followup_barcode'):
                    follow_up_barcode = c_text
                follow_up_data = { 'sequence' : sequence}
                for c_node, c_stack, c_attr, c_text in dom_scan(s_node, 'follow_up/*'):
                    if 'xsd_ver' in c_attr:
                        follow_up_data[c_attr.get('preferred_name', c_stack[-1])] = c_text
                self.emit( state, follow_up_barcode, follow_up_data, 'Followup' )

        return state

    def emit(self, state, key, entry, entry_type):
        state['generators'][entry_type].find(entry)

def initial_state():
    individual_generator = IndividualGenerator()
    biosample_generator = BiosampleGenerator(individual_generator.gid)

    return record_initial_state({'Individual': individual_generator, 'Biosample': biosample_generator})

def build_processor(extract, subtype):
    def process(state, file):
        raw = file.read()
        try:
            print('parsing', file.name)
            dom = parseString(raw)
        except:
            dom = None
                        
        if dom:
            return extract.parseXMLFile(state, dom, subtype)
        else:
            return state

    return process

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('subtype')
    # parser.add_argument('--clinical', type=str, help='path to directory containing clinical files')
    # parser.add_argument('--biospecimen', type=str, help='path to ')
    parser.add_argument('--input', type=str, help='path to directory containing the input files')
    parser.add_argument('--output', type=str, help='path to output file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    extract = ClinicalParser()
    process = build_processor(extract, args.subtype)
    
    state = initial_state()
    process_input(state, args.input, process)
    output_state(state, args.output)

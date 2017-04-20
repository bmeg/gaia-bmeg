#!/usr/bin/env python

__all__ = ['DrugBankReader']

import sys
import re
import xml.sax
import json
import logging

reWord = re.compile(r'\w')
reSpace = re.compile(r'\s')


def ignore(e, v,**kwds):
    return None

def create_list(e, v):
    return [v]

def string_pass(e, v):
    return v

def pass_data(e,v,**kwds):
    return kwds

def create_dict_list(e, v, **kwds):
    return [kwds]

def debug_emit(e, v, **kwds):
    print json.dumps(kwds)
    

f_map = [
    (['drugbank','drug'], None, debug_emit),    

    (['drugbank','drug','drugbank-id'], None, create_list),
    
    (['drugbank','drug',
        ['name','description','cas-number','unii','state','indication', 'toxicity','metabolism','half-life',
        'clearance', 'absorption','mechanism-of-action','route-of-elimination','average-mass','monoisotopic-mass','protein-binding']], 
        None, string_pass),

    (['drugbank','drug','synonyms','synonym'], None, create_list),
    (['drugbank','drug','synonyms'], None, pass_data),

    (['drugbank','drug','products','product','*'], None, string_pass),
    (['drugbank','drug','products','product'], 'products', create_dict_list),
    (['drugbank','drug','products'], None, pass_data),
    
    (['drugbank','drug','groups', 'group'], None, create_list),
    (['drugbank','drug','groups'], None, pass_data),
    

    (['drugbank','drug','targets','target','polypeptide','go-classifiers','go-classifier', '*'], None, string_pass),    
    (['drugbank','drug','targets','target','polypeptide','go-classifiers','go-classifier'], 'go-classifiers', create_dict_list),
    (['drugbank','drug','targets','target','polypeptide','go-classifiers'], None, pass_data),
    
    (['drugbank','drug','targets','target','polypeptide','synonyms','synonym'], None, create_list),
    (['drugbank','drug','targets','target','polypeptide','synonyms'], None, pass_data),
    
    (['drugbank','drug','targets','target','polypeptide',
        ['name','general-function','specific-function','gene-name','locus','cellular-location','molecular-weight','organism']], None, pass_data),
    
    (['drugbank','drug','targets','target','polypeptide'], None, pass_data),

    (['drugbank','drug','targets','target','references','articles','article','*'], None, string_pass),
    (['drugbank','drug','targets','target','references','articles','article'], 'articles', create_dict_list),
    (['drugbank','drug','targets','target','references','articles'], None, pass_data),
    (['drugbank','drug','targets','target','references'], None, pass_data),
    
    (['drugbank','drug','enzymes','enzyme','polypeptide','go-classifiers','go-classifier','*'], None, string_pass ),
    (['drugbank','drug','enzymes','enzyme','polypeptide','go-classifiers','go-classifier'], 'go-classifiers', create_dict_list ),
    (['drugbank','drug','enzymes','enzyme','polypeptide','go-classifiers'], None, pass_data ),
    
    (['drugbank','drug','enzymes','enzyme','polypeptide',
        ['name','general-function','specific-function','gene-name','organism','molecular-weight','amino-acid-sequence','gene-sequence']], None, string_pass ),
    
    (['drugbank','drug','enzymes','enzyme','polypeptide'], None, pass_data ),
    (['drugbank','drug','enzymes','enzyme'], 'enzymes', create_dict_list ),
    (['drugbank','drug','enzymes'], None, pass_data ),

    (['drugbank','drug','targets','target',
        ['id','name','organism','known-action']], None, pass_data),

    (['drugbank','drug','targets','target'], 'targets', create_dict_list),
    (['drugbank','drug','targets'], None, pass_data),

    (['drugbank','drug','classification','description', '*'], None, string_pass),
    (['drugbank','drug','classification','description'], None, pass_data),
    (['drugbank','drug','classification',['direct-parent','kingdom','superclass','class','subclass']], None, string_pass),
    (['drugbank','drug','classification'], None, pass_data),
    
    (['drugbank','drug','categories','category', '*'], None, string_pass),
    (['drugbank','drug','categories','category'], 'categories', create_dict_list),
    (['drugbank','drug','categories'], None, pass_data),

    (['drugbank','drug','mixtures','mixture','*'], None, string_pass),
    (['drugbank','drug','mixtures','mixture'], 'mixtures', create_dict_list),
    (['drugbank','drug','mixtures'], None, pass_data),
    
    (['drugbank','drug','drug-interactions','drug-interaction','*'], None, string_pass),
    (['drugbank','drug','drug-interactions','drug-interaction'], 'drug-interactions', create_dict_list),
    (['drugbank','drug','drug-interactions'], None, pass_data),

    (['drugbank','drug','dosages','dosage','*'], None, string_pass),
    (['drugbank','drug','dosages','dosage'], 'dosages', create_dict_list),
    (['drugbank','drug','dosages'], None, pass_data),


]

class StackLevel:
    def __init__(self, name):
        self.data = {}
        self.name = name
        self.has_children = False

def stack_match(query, elem):
    match = True
    if len(query) != len(elem):
        return False
    for q, e in zip(query, elem):
        if isinstance(q, list):
            if e not in q:
                return False
        else:
            if q != "*" and q != e:
                return False
    return True

NOT_FOUND = set()

class DrugBankHandler(xml.sax.ContentHandler):
    def __init__(self, record_write):
        xml.sax.ContentHandler.__init__(self)
        self.record_write = record_write
        self.stack = []

    def startElement(self, name, attrs):
        if len(self.stack):
            self.stack[-1].has_children = True
        self.stack.append(StackLevel(name))
        self.buffer = ""

    def characters(self, text):
        self.buffer += text

    def endElement(self, name):
        found = False
        stack_id = list(i.name for i in self.stack)
        #print stack_id
        level = self.stack.pop()
        for s, out_name, f in f_map:
            if stack_match(s, stack_id):
                found = True
                if out_name is None:
                    out_name = stack_id[-1]
                v = f(self.record_write, self.buffer, **level.data)
                if v is not None:
                    if isinstance(v,list):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name].extend(v)
                        else:
                            self.stack[-1].data[out_name] = v
                    elif isinstance(v,dict):
                        if out_name in self.stack[-1].data:
                            self.stack[-1].data[out_name] = dict(self.stack[-1].data, **v)
                        else:
                            self.stack[-1].data[out_name] = v
                    else:
                        self.stack[-1].data[out_name] = v
        if not found:
            n = ",".join(stack_id)
            if n not in NOT_FOUND:
                NOT_FOUND.add(n)
                logging.warning("combiner for %s not found" % (",".join(stack_id)))
        self.buffer = ""
        


def run_sax_parse(handle, queue):
    def record_write(record):
        queue.send(record)
        
    handler = DrugBankHandler(record_write)
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    parser.parse(handle)
    queue.close()

def emit(msg):
    print msg

def parse_drugbank(handle):   
    handler = DrugBankHandler(emit)
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    parser.parse(handle)

if __name__ =="__main__":
    logging.basicConfig(level=logging.INFO)
    path = sys.argv[1]
    handle = open(path)

    parse_drugbank(handle)
    
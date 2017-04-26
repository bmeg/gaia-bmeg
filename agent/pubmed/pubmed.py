#!/usr/bin/env python

import sys
import re
import xml.sax
import json
import logging
import gzip 

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

class PubMedHandler(xml.sax.ContentHandler):
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
        
def emit(msg):
    print msg

def parse_pubmed(handle):   
    handler = PubMedHandler(emit)
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    parser.parse(handle)

if __name__ =="__main__":
    logging.basicConfig(level=logging.INFO)
    path = sys.argv[1]
    handle = gzip.open(path)

    parse_pubmed(handle)
    
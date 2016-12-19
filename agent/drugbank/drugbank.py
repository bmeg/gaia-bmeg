#!/usr/bin/env python

__all__ = ['DrugBankReader']


import re
import xml.sax
import json
import multiprocessing


reWord = re.compile(r'\w')
reSpace = re.compile(r'\s')

ont_base = "http://bmeg.ucsc.edu/drug_ont#"

def id_mutate(s):
    return "http://bmeg.ucsc.edu/drugbank/" + s

def default_mutate(s):
    return  {"@value" : s }

def link_mutate(s):
    return  {"@id" : s }

def id_link_mutate(s):
    return link_mutate(id_mutate(s))

def ont_mutate(s):
    return link_mutate(ont_base + reSpace.sub("_", s))

def m_copy(s):
    return s

extenal_mapping = {
    "KEGG Drug" : "http://bmeg.ucsc.edu/kegg/",
    "KEGG Compound" : "http://bmeg.ucsc.edu/kegg/",
    "Drugs Product Database (DPD)" : "http://bmeg.ucsc.edu/DPD/",
    "National Drug Code Directory" : "http://bmeg.ucsc.edu/NCDC/",
    "PharmGKB" : "http://bmeg.ucsc.edu/PharmGKB/",
    "UniProtKB" : "http://bmeg.ucsc.edu/uniprotkb/",
    "GenBank" : "http://bmeg.ucsc.edu/genebank/",
    "ChEBI" : "http://bmeg.ucsc.edu/chebi/",
    "BindingDB" : "http://bmeg.ucsc.edu/BindingDB/",
    "IUPHAR" : "http://bmeg.ucsc.edu/IUPHAR/",
    "Guide to Pharmacology" : "http://bmeg.ucsc.edu/g2pharma/",
    "PubChem Compound" : "http://bmeg.ucsc.edu/pubchem/",
    "PubChem Substance" : "http://bmeg.ucsc.edu/pubchem/",
    "ChemSpider" : "http://bmeg.ucsc.edu/ChemSpider/",
    "PDB" : "http://bmeg.ucsc.edu/pdb/",
}

def external_ident(s):
    return extenal_mapping[s['resource']] + s['identifier']

def property_select(s):
    return ont_base + reSpace.sub("_", s['kind'])

def property_format(s):
    return s['value']

f_map = [
    (['drugs', 'drug', 'drugbank-id'], "@id", id_mutate, False),
    (['drugs', 'drug', 'cas-number'],  ont_base + "casNumber", None, True),

    (['drugs', 'drug', 'half-life'], ont_base + "halfLife", None, True),

    (['drugs', 'drug', 'transporters', 'transporter', 'actions', 'action'], ont_base + "transporterAction", ont_mutate, True),
    (['drugs', 'drug', 'targets', 'target', 'actions', 'action'], ont_base + "targetAction", ont_mutate, True),
    (['drugs', 'drug', 'groups', 'group'], ont_base + "group", None, True),
    (['drugs', 'drug', 'enzymes', 'enzyme', 'actions', 'action'], ont_base + "group", ont_mutate, True),
    (['drugs', 'drug', 'clearance'], ont_base + "clearance", None, True),
    (['drugs', 'drug', 'categories', 'category'], ont_base + "category", ont_mutate, True),
    (['drugs', 'drug', 'synonyms', 'synonym'], ont_base + "synonym", None, True),
    (['drugs', 'drug', 'name'], ont_base + "name", None, True),
    (['drugs', 'drug', 'brands', 'brand'], ont_base + "brand", None, True),
    (['drugs', 'drug', 'targets', 'target', 'references'], ont_base + "targetReference", None, True),

    (['drugs', 'drug', 'description'], ont_base + "description", None, True),

    (['drugs', 'drug', 'general-references'],  ont_base + "generalReference", None, True),
    (['drugs', 'drug', 'external-links', 'external-link', 'url'], ont_base + "xref", link_mutate, True),

    (['drugs', 'drug', 'indication'], ont_base + "indication", None, True),

    (['drugs', 'drug', 'calculated-properties', 'property', 'kind'], "kind", m_copy, False),
    (['drugs', 'drug', 'calculated-properties', 'property', 'value'], "value", m_copy, False),
    (['drugs', 'drug', 'calculated-properties', 'property'], property_select, property_format, True),


    (['drugs', 'drug', 'drug-interactions', 'drug-interaction', 'drug'], ont_base + "drugInteraction", id_link_mutate, True),
    (['drugs', 'drug', 'drug-interactions', 'drug-interaction', 'name'], ont_base + "name", None, True),
    (['drugs', 'drug', 'drug-interactions', 'drug-interaction', 'description'], ont_base + "description", None, True),
    (['drugs', 'drug', 'drug-interactions', 'drug-interaction'], ont_base + "drugInteraction", m_copy, True),

    (['drugs', 'drug', 'external-identifiers', 'external-identifier', 'resource'], "resource", m_copy, False),
    (['drugs', 'drug', 'external-identifiers', 'external-identifier', 'identifier'], "identifier", m_copy, False),
    (['drugs', 'drug', 'external-identifiers', 'external-identifier'], ont_base + "extIdent", external_ident, True)
]

class StackLevel:
    def __init__(self, name):
        self.data = {}
        self.name = name
        self.has_children = False

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
        if not level.has_children:
            for s, f, m, append in f_map:
                if s == stack_id:
                    found = True
                    if f not in self.stack[-1].data:
                        self.stack[-1].data[f] = []

                    val = None
                    if m is None:
                        val = default_mutate(self.buffer)
                    else:
                        val = m(self.buffer)

                    if append:
                        if len(val):
                            self.stack[-1].data[f].append( val )
                    else:
                        if len(val):
                            self.stack[-1].data[f] = val
            #if not found:
            #    sys.stderr.write("Not Found: %s %s\n" % (str(stack_id), self.buffer.encode('ascii', errors='ignore')))
        else:
            found = False
            for s, f, m, append in f_map:
                if s == stack_id:   
                    found = True
                    if callable(f):
                        f = f(level.data)
                    if f not in self.stack[-1].data:
                        self.stack[-1].data[f] = []

                    val = m(level.data)
                    if append:
                        if len(val):
                            self.stack[-1].data[f].append( val )
                    else:
                        if len(val):
                            self.stack[-1].data[f] = val

            if not found:
                if len(self.stack):
                    for k,v in level.data.iteritems():
                        self.stack[-1].data[k] = v
        
            
        self.buffer = ""
        post = list(i.name for i in self.stack)
        #print "post", post
        if post == ['drugs']:
            #print self.stack[-1].data
            self.record_write(self.stack[-1].data)
            self.cur_record = {}


def run_sax_parse(handle, queue):
    def record_write(record):
        queue.send(record)
        
    handler = DrugBankHandler(record_write)
    parser = xml.sax.make_parser()
    #parser.setFeature(xml.sax.handler.feature_namespaces, True)
    parser.setContentHandler(handler)
    parser.parse(handle)
    queue.close()


def parse_drugbank(handle):   
    q_out, q_in = multiprocessing.Pipe(duplex = False)
    process = multiprocessing.Process(target=run_sax_parse, args=(handle, q_in))
    process.start()
    q_in.close()

    while True:
        try:
            yield q_out.recv()
        except EOFError:
            break


if __name__ =="__main__":
    import sys
    path = sys.argv[1]
    handle = open(path)

    for line in parse_drugbank(handle):
        #pass 
        print json.dumps(line, indent=4)
    
#!/usr/bin/env python

# parses GO annotation file
# http://www.geneontology.org/doc/GO.references

import sys
import csv
import re
import json

"""
bmeg_ns = Namespace("http://purl.org/bmeg/owl#")
hugo_ns = Namespace("http://genenames.org/genes/")
genefamily_ns = Namespace("http://genenames.org/genefamilies/")
locus_ns = Namespace("http://purl.org/bmeg/gene_locus/")
tcga_samples = Namespace("http://purl.org/bmeg/tcga/")
pubmed_ns = Namespace("http://identifiers.org/pubmed/")
gene_ns = Namespace("http://identifiers.org/insdc/")
entrez_ns = Namespace("http://identifiers.org/ncbigene/")
mgd_ns = Namespace("http://identifiers.org/mgd/")
rgd_ns = Namespace("http://identifiers.org/rgd/")
go_ns = Namespace("http://purl.org/obo/owl/GO#")
uniprot_ns = Namespace("http://identifiers.org/uniprot/")
obo = Namespace("http://www.geneontology.org/formats/oboInOwl#")

goref_ns = Namespace("http://purl.org/obo/owl/GO_REF#")
"""

def uniprot_ns(n):
	return "uniprot:" + n

def go_ns(n):
	return "go:" + n

def goref_ns(n):
	return "goref:" + n

def gene_ns(n):
	return "gene:" + n

def pubmed_ns(n):
	return "pubmed:" + n

UNIPROT_COL = 1
SYMBOL_COL = 2
GOID_COL = 4
REF_COL = 5
EVIDENCE_COL = 6
NAME_COL = 9

if __name__ == "__main__":

	for line in sys.stdin:
		if not line.startswith("!"):
			row = line.rstrip().rsplit("\t")
			go_id = row[GOID_COL].replace("GO:", "")
			node = {}
			node['class'] = "GeneFunctionAnnotation"
			node['geneEdges'] = []
			for r in row[SYMBOL_COL].split("|"):
				if re.search(r'^[A-Z0-9]+$', r):
					node['geneEdges'].append( gene_ns(r) ) 
			node['geneEdges'].append( uniprot_ns(row[UNIPROT_COL]) )
			node['functionEdges'] = [ go_ns(go_id) ]
			node['evidenceEdges'] = [ go_ns(row[EVIDENCE_COL]) ]

			node['referenceEdges'] = []
			if row[REF_COL].startswith("GO_REF:"):
				node['referenceEdges'].append( goref_ns(row[REF_COL].replace("GO_REF:", "")) )
			if row[REF_COL].startswith("PMID:"):
				node['referenceEdges'].append( pubmed_ns(row[REF_COL].replace("PMID:", "")))

			if len(row[NAME_COL]):
				node['title'] = row[NAME_COL] 
			
			print json.dumps(node)

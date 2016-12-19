#!/usr/bin/env python

import sys
import json
import argparse
import urllib
from bmeg import phenotype_pb2
from google.protobuf import json_format

BASE_URL = "http://lincs.hms.harvard.edu/db/api/v1/"

def action_parse(args):
    with open(args.file) as handle:
        data = json.load(handle)
    
    dataRecords = {}
    for a in data:
        key = a['datarecordID']
        if key not in dataRecords:
            dataRecords[key] = {}
        for b in ["smallmolecule_smName", 'smallmolecule_smLincsID', 'cell_clName' ]:
            dataRecords[key][b] = a[b]
        dataRecords[key][a['datapointName']] = a['datapointValue']
        
    cell_line_vs_drug = {}
    for v in dataRecords.values():
        c_name = v['cell_clName']
        sm_name = v['smallmolecule_smName']
        if c_name not in cell_line_vs_drug:
            cell_line_vs_drug[c_name] = {}
        if sm_name not in cell_line_vs_drug[c_name]:
            cell_line_vs_drug[c_name][sm_name] = []
        """
        cell_line_vs_drug[c_name][sm_name].append({
            'smallMolConcentration' : v['smallMolConcentration'],
            'relativeCellCount' : v['relativeCellCount'],
            'timepoint' : v['timepoint'],
            'relativeGrowthRateAfterTreatment' : v['relativeGrowthRateAfterTreatment'],
            'totalControlCellCount' : v['totalControlCellCount'],
            'replicate' : v['replicate']
        })
        """
        cell_line_vs_drug[c_name][sm_name].append(v)

    #for k, v in dataRecords.items():
    #    print k, json.dumps(v)
    out = []
    for cell in cell_line_vs_drug:
        for drug in cell_line_vs_drug[cell]:
            rc = phenotype_pb2.ResponseCurve()
            rc.sample = cell
            rc.compound = drug
            rc.responseType = rc.GROWTH
            for i in cell_line_vs_drug[cell][drug]:
                t = rc.values.add()
                t.dose = float(i['smallMolConcentration'])
                t.response = float(i['relativeCellCount'])
            out.append(rc)
    
    for i in out:
        print json_format.MessageToJson(i)
    #print cell_line_vs_drug

def action_list(args):
    
    handle = urllib.urlopen(BASE_URL + "dataset")
    print handle.read()
    handle.close()
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    subparsers = parser.add_subparsers(title="subcommand")
    
    parser_list = subparsers.add_parser('list')
    parser_list.set_defaults(func=action_list)
    
    parser_parse = subparsers.add_parser("parse")
    parser_parse.add_argument("file")
    parser_parse.set_defaults(func=action_parse)

    args = parser.parse_args()
    func = args.func
    e = func(args)
    sys.exit(e)
    
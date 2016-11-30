#!/usr/bin/env python

import sys
import json


def imaging_import(file):
    with open(file) as handle:
        data = json.load(handle)
    
    dataRecords = {}
    for a in data:
        key = a['datarecordID']
        if key not in dataRecords:
            dataRecords[key] = {}
        for b in ["smallmolecule_smName", 'smallmolecule_smLincsID', 'cell_clName' ]:
            dataRecords[key][b] = a[b]
        dataRecords[key][a['datapointName']] = a['datapointValue']

    for k, v in dataRecords.items():
        print k, json.dumps(v)

if __name__ == "__main__":
    imaging_import(sys.argv[1])
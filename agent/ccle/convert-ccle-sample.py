#!/usr/bin/env python

import json
import sys
import csv

"""

Clinical data download
curl -o CCLE_sample_info_file_2012-10-18.txt "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801"

"""

if __name__ == "__main__":

    with open(sys.argv[1]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line in reader:
            
            out = {
                "type" : "Sample",
                "gid" : "ccle:%s" % (line['CCLE name']),
                "sex" : line['Gender'],
                "histology" : line["Histology"],
                "tissue" : line["Site Primary"],
                "subtype" : line["Hist Subtype1"],
                "alias" : line["Cell line primary name"],
                "source" : line["Source"],
                "notes" : line["Notes"]
            }
            #print json.dumps(line)
            print json.dumps(out)
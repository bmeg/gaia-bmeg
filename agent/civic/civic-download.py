#!/usr/bin/env python

import json
import requests


def get_genes():
    page = 1
    while True:
        data = requests.get("https://civic.genome.wustl.edu/api/genes", params={"page" : page}).json()
        for r in data['records']:
            yield r
        if data["_meta"]["current_page"] >= data["_meta"]["total_pages"]:
            break
        page += 1

def get_variant(variant_id):
    data = requests.get("https://civic.genome.wustl.edu/api/variants/%s" % (variant_id)).json()
    return data

found_variants = set()

with open("civic-genes.out", "w") as handle:
    for g in get_genes():
        handle.write("%s\n" % json.dumps(g))
        for i in g['variants']:
            found_variants.add(i['id'])

with open("civic-variants.out", "w") as handle:
    for i in found_variants:
        o = get_variant(i)
        handle.write("%s\n" % json.dumps(o))
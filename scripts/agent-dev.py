#!/usr/bin/env python

import os
import argparse
from glob import glob
import yaml
import hashlib
import urllib
import tempfile
import subprocess
import json

CACHE_DIR=os.path.join( os.path.dirname( os.path.abspath(__file__)), "..", "cache" )
WORK_DIR=os.path.join( os.path.dirname( os.path.abspath(__file__)), "..", "work" )

def download_cache(url):
    u =  hashlib.md5(url).hexdigest()
    f = os.path.join(CACHE_DIR, u)
    if not os.path.exists(f):
        urllib.urlretrieve(url, f)
    return f

def which(file):
    for path in os.environ["PATH"].split(":"):
        p = os.path.join(path, file)
        if os.path.exists(p):
            return p

def prep_inputs(config):
    inputs = {}
    for k,v in config["inputs"].items():
        if v['type'] == "Download":
            url = v['sourceUrl']
            path = download_cache(url)
            inputs[k] = {
                "class" : "File",
                "path" : path
            }
    return inputs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("dir")
    
    args = parser.parse_args()
    
    for agent in glob(os.path.join(args.dir, "*", "*.agent.yml")):
        with open(agent) as handle:
            config = yaml.load(handle)
            print config
            inputs = prep_inputs(config)
            workflow = os.path.join(os.path.dirname(os.path.abspath(agent)), config['workflow']['run'] )
            print workflow
            
            workflow_inputs = {}
            for k, v in config['workflow']['inputs'].items():
                workflow_inputs[k] = inputs[v]
            
            out = tempfile.NamedTemporaryFile(dir=WORK_DIR, prefix="agent_workflow_", suffix=".json", delete=False)
            out_name = out.name
            out.write(json.dumps(workflow_inputs))
            out.close()

            cmd = [
                which("cwltool"),
                "--tmp-outdir-prefix=" + WORK_DIR + "/tmp",
                "--tmpdir-prefix=" + WORK_DIR + "/tmp",
                workflow,
                out_name
            ]

            print "Running %s" % (" ".join(cmd))
            print("workdir: " + WORK_DIR)
            proc = subprocess.Popen(cmd, cwd=WORK_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            print "stdout:", stdout
            print "stderr:", stderr
            outputs = json.loads(stdout)
            
            for o in config['workflow']['outputs']:
                print o, outputs[o]
            
    

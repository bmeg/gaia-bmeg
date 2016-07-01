#!/usr/bin/env bash

# Sample ctdd converter script. Please set-up with your local paths and options.

CWD=$(pwd)
CONVERTER_SCRIPT=${CWD}/convert-ctdd.py

OUTPATH="${CWD}/tmp_converted/ctdd_converted.json"

RESPONSE="${CWD}/raw-data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt"
METADRUG="${CWD}/raw-data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt"
METACCL="${CWD}/raw-data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt"
METAEXPERIMENT="${CWD}/raw-data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt"


python $CONVERTER_SCRIPT \
--response $RESPONSE \
--metadrug $METADRUG \
--metacellline $METACCL \
--metaexperiment $METAEXPERIMENT \
--out $OUTPATH \
--format json

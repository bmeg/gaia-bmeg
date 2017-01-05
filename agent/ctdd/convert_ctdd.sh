#!/usr/bin/env bash

# Sample ctdd converter script. Please set-up with your local paths and options.

CWD=$(pwd)
CONVERTER_SCRIPT=${CWD}/convert-ctdd.py

DATADIR=$1

OUTPATH="${CWD}/ctdd_converted.json"

RESPONSE="${DATADIR}/v20.data.curves_post_qc.txt"
METADRUG="${DATADIR}/v20.meta.per_compound.txt"
METACCL="${DATADIR}/v20.meta.per_cell_line.txt"
METAEXPERIMENT="${DATADIR}/v20.meta.per_experiment.txt"
AVEDATA="${DATADIR}/v20.data.per_cpd_avg.txt"

python $CONVERTER_SCRIPT \
--response $RESPONSE \
--metadrug $METADRUG \
--metacellline $METACCL \
--metaexperiment $METAEXPERIMENT \
--data $AVEDATA \
--multi $OUTPATH \
--format json

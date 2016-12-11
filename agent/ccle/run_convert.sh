#!/bin/bash

BDIR=$( cd $(dirname $0) && pwd )

CONVERTER_SCRIPT=${BDIR}/convert-ccle.py

OUTPATH="ccle_converted.json"

MAF=$1/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf
CSV=$1/CCLE_NP24.2009_Drug_data_2015.02.24.csv
EXP=$1/CCLE_Expression_2012-09-29.res
SAM=$1/CCLE_sample_info_file_2012-10-18.txt

python $CONVERTER_SCRIPT --maf $MAF --out ccle.mutations.json --format json
python $CONVERTER_SCRIPT --drug $CSV --out ccle.drug.json --format json
python $CONVERTER_SCRIPT --sample $SAM --out ccle.sample.json --format json
python $CONVERTER_SCRIPT --exp $EXP --out ccle.expression.json --format json

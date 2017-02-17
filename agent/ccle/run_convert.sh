#!/usr/bin/env bash

#To Run
#curl -o CCLE_NP24.2009_Drug_data_2015.02.24.csv "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_Drug_data_2015.02.24.csv?downloadff=true&fileId=20777"
#curl -o CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_26/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz?downloadff=true&fileId=6873"
#curl -o CCLE_sample_info_file_2012-10-18.txt "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801"
#curl -o CCLE_Expression_2012-09-29.res "https://portals.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_21/CCLE_Expression_2012-09-29.res?downloadff=true&fileId=6760"
#./run_convert.sh ./

BDIR=$( cd $(dirname $0) && pwd )

CONVERTER_SCRIPT=${BDIR}/convert-ccle.py


MAF=$1/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf
CSV=$1/CCLE_NP24.2009_Drug_data_2015.02.24.csv
EXP=$1/CCLE_Expression_2012-09-29.res
SAM=$1/CCLE_sample_info_file_2012-10-18.txt

python $CONVERTER_SCRIPT --drug $CSV --out ccle.ResponseCurve.json --format json
python $CONVERTER_SCRIPT --sample $SAM --out ccle.Biosample.json --format json
python $CONVERTER_SCRIPT --exp $EXP --out ccle.GeneExpression.json --format json

#!/bin/bash

# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17_fitted_dose_response.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/GDSC-CCLE-CTRP_conversion.xlsx
# curl -O ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Screened_Compounds.xlsx

./gdsc-convert.py GDSC-CCLE-CTRP_conversion.xlsx Cell_Lines_Details.xlsx Screened_Compounds.xlsx v17_fitted_dose_response.xlsx


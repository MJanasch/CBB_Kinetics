#!/usr/bin/env bash

# Pipeline for extracting and plotting SKM metabolite concentrations

source $1

### 1. EXTRACT METABOLITES HEADER FROM MATLAB FILE #############################
matlab -nojvm -r "infile='$MODEL_FILE'" < $HEADER_SCRIPT > /dev/null
HEADER=$(cat /tmp/cbb_metabolite_header.long.txt | grep -v "^BioMass_" | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

### 2. EXTRACT CONCENTRATIONS FROM MATLAB FILES ################################
matlab -nojvm -r "indir='$CONC_DIR'" < $CONC_SCRIPT > /dev/null

### 3. CONSTRUCT CONCENTRATIONS FILE ###########################################
CONC_FILE="${CONC_DIR}/all_metabolite_concentrations.tab"
(echo -e "$HEADER"; cat /tmp/cbb_metabolite_concs.tab) > $CONC_FILE

### 4. PLOT CONCENTRATIONS #####################################################
$PLOT_SCRIPT $CONC_FILE

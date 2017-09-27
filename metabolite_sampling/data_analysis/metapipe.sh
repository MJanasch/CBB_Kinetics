#!/usr/bin/env bash

# Pipeline for extracting and plotting SKM metabolite concentrations

source $1

### 1. EXTRACT METABOLITES HEADER FROM MATLAB FILE #############################

# Report progress
echo -e "\n\e[94mStep 1: Extracting metabolites header...\e[0m\n"

matlab -nojvm -r "infile='$MODEL_FILE'" < $HEADER_SCRIPT > /dev/null
HEADER=$(cat /tmp/cbb_metabolite_header.long.txt | grep -v "^BioMass_" | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Report step done
echo -e "\n\e[92mStep 1: Done.\e[0m\n"

### 2. EXTRACT CONCENTRATIONS FROM MATLAB FILES ################################

# Report progress
echo -e "\n\e[94mStep 2: Extracting concentrations...\e[0m\n"

matlab -nojvm -r "indir='$CONC_DIR'" < $CONC_SCRIPT > /dev/null

# Report step done
echo -e "\n\e[92mStep 2: Done.\e[0m\n"

### 3. CONSTRUCT CONCENTRATIONS FILE ###########################################

# Report progress
echo -e "\n\e[94mStep 3: Finalizing concentration file...\e[0m\n"

CONC_FILE="${CONC_DIR}/all_metabolite_concentrations.tab"
(echo -e "$HEADER"; cat /tmp/cbb_metabolite_concs.tab) > $CONC_FILE

rm /tmp/cbb_metabolite_header.long.txt /tmp/cbb_metabolite_concs.tab

# Report step done
echo -e "\n\e[92mStep 3: Done.\e[0m\n"

### 4. PLOT CONCENTRATIONS #####################################################

# Report progress
echo -e "\n\e[94mStep 4: Plotting concentrations...\e[0m\n"

$PLOT_SCRIPT $CONC_FILE

# Report step done
echo -e "\n\e[92mStep 4: Done.\e[0m\n"

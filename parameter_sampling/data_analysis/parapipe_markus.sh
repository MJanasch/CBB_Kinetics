#!/usr/bin/env bash

# Pipeline for extracting and plotting SKM parameter sets

source $1

### 1. CREATE LIST WITH PERCENT STABLE STEADY STATES PER CONCENTRATION SET #####

# Report progress
echo -e "\n\e[94mStep 1: Calculating percent stable steady states...\e[0m\n"

# Calculate percent stable steady states
matlab -nojvm -r "indir='$PARA_DIR'" < $STAB_SCRIPT > /dev/null

### Output here is 'met_set_vs_percent_steady.tab', a table with one column of
### metabolite set number and a second column of % stable steady states

# Report step done
echo -e "\n\e[92mStep 1: Done.\e[0m\n"

### 2. EXTRACT FCC VALUES ######################################################

# Report progress
echo -e "\n\e[94mStep 2: Extracting flux control coefficients...\e[0m\n"

# Clean up temporary output directory
if [ -d $TMP_DIR ]; then
  rm -r $TMP_DIR
fi

# Create temporary output directory
mkdir -p $TMP_DIR/fcc

# Extract FCCs
m_vars="indir='$PARA_DIR'; model_file='$MODEL_FILE'; outdir='$TMP_DIR/fcc/'"
matlab -nojvm -r "$m_vars" < $FCCEX_SCRIPT > /dev/null

### Output here is a tab file for every stable steady state for every metabolite set, containing
### a column with consisting only of metabolite set number, second column is only the number of 
### the current stable state, third column is reaction numbers from 1 to the end, then comes the 
### actual FCC vector

# Create FCCs header
matlab -nojvm -r "model_file='$MODEL_FILE'; outdir='$TMP_DIR/';" \
< $HEADER_SCRIPT > /dev/null

### Output here is a .txt file containing all the names of the reactions (as a column)


HEADER=$(echo -en "Conc_set\tStable_set\tReaction\t"; \
cat $TMP_DIR/cbb_reaction_header.long.txt | grep -v "^BioMass_" | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Concatenate FCC data
(echo -e "$HEADER"; find $TMP_DIR/fcc/ -name "*.tab" | xargs -n 32 -P 8 cat) > \
$PARA_DIR/concset_stabstate_rxn_FCCs.tab

### Output here is a .tab file containing all FCC data, first column is the met set number, the
### second column is the number of the stable steady state of the particular met set, the third
### column is the reaction number (1 to the end, 27 or 30), the other columns are the reactions
### and their associated FCC (Header is Conc_set, Stable_set, Reaction, RuBisCO, PGK,...)


# Gzip and retire
pigz $PARA_DIR/concset_stabstate_rxn_FCCs.tab
retire $PARA_DIR/concset_stabstate_rxn_FCCs.tab.gz

# Report step done
echo -e "\n\e[92mStep 2: Done.\e[0m\n"


### 3. EXTRACT CCC VALUES ######################################################

# Report progress
echo -e "\n\e[94mStep 3: Extracting concentration control coefficients...\e[0m\n"

# Create temporary output directory
mkdir -p $TMP_DIR/ccc

# Extract CCCs
m_vars="indir='$PARA_DIR'; model_file='$MODEL_FILE'; outdir='$TMP_DIR/ccc/'"
matlab -nojvm -r "$m_vars" < $CCCEX_SCRIPT > /dev/null

# Create CCCs header
matlab -nojvm -r "model_file='$MODEL_FILE'; outdir='$TMP_DIR/';" \
< $MET_NAMES_SCRIPT > /dev/null

### Output here is a .txt file containing all the names of the reactions (as a column)


MET_HEADER=$(echo -en "Conc_set\tStable_set\tMetabolite\t"; \
cat $TMP_DIR/cbb_reaction_header.long.txt | grep -v "^BioMass_" | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Concatenate CCC data
(echo -e "$MET_HEADER"; find $TMP_DIR/ccc/ -name "*.tab" | xargs -n 32 -P 8 cat) > \
$PARA_DIR/concset_stabstate_met_CCCs.tab

# Gzip and retire
pigz $PARA_DIR/concset_stabstate_met_CCCs.tab
retire $PARA_DIR/concset_stabstate_met_CCCs.tab.gz

# Report step done
echo -e "\n\e[92mStep 3: Done.\e[0m\n"



### 4. EXTRACT SAMPLED PARAMETERS ##############################################

# Report progress
echo -e "\n\e[94mStep 4: Extracting parameters...\e[0m\n"

# Create temporary output directory
mkdir -p $TMP_DIR/par

# Extract parameters
m_vars="indir='$PARA_DIR'; outdir='$TMP_DIR/par/'"
matlab -nojvm -r "$m_vars" < $PAREX_SCRIPT > /dev/null

# Store parameters header (bulk created by PAREX_SCRIPT above)
HEADER=$(echo -en "Conc_set\tStable\t"; cat $TMP_DIR/par/par_header.long.txt | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Concatenate parameters data
(echo -e "$HEADER"
find $TMP_DIR/par/ -name "*.tab" | while read File; do
  Conc_set=`echo $File | rev | cut -f 1 -d \/ | rev | cut -f 1 -d \.`
  cat $File | sed -e "s/^/${Conc_set}\t/"
done) > $PARA_DIR/concset_stability_parameters.tab

# Gzip and retire
pigz $PARA_DIR/concset_stability_parameters.tab
retire $PARA_DIR/concset_stability_parameters.tab.gz

# Report step done
echo -e "\n\e[92mStep 4: Done.\e[0m\n"

### 4. PLOT METABOLITE CONCENTRATIONS AND RATIOS VS STABILITY ##################

# Plot with R
$CSTAB_SCRIPT $CONC_FILE $PARA_DIR/met_set_vs_percent_steady.tab

### 5. PLOT FCC HEATMAP ########################################################

# Plot with R
$HTMAP_SCRIPT $PARA_DIR/concset_stabstate_rxn_FCCs.tab.gz \
$TMP_DIR/cbb_reaction_header.long.txt

$HTMAP_CCC_SCRIPT $PARA_DIR/concset_stabstate_met_CCCs.tab.gz \
$TMP_DIR/cbb_met_list.long.txt

### 6. PLOT KM VS CONCENTRATION ################################################

### 7. DIMENSIONALITY REDUCTION OF KM OVER CONCENTRATION #######################

### 8. CLUSTERING OF KM OVER CONCENTRATION #####################################

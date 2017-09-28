#!/usr/bin/env bash

# Pipeline for extracting and plotting SKM parameter sets

source $1

### 1. CREATE LIST WITH PERCENT STABLE STEADY STATES PER CONCENTRATION SET #####

# Report progress
echo -e "\n\e[94mStep 1: Calculating percent stable steady states...\e[0m\n"

# Calculate percent stable steady states
matlab -nojvm -r "indir='$PARA_DIR'" < $STAB_SCRIPT > /dev/null

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
m_vars="indir='$PARA_DIR'; model_file='$MODEL_FILE'; outdir='/tmp/skm/fcc/'"
matlab -nojvm -r "$m_vars" < $FCCEX_SCRIPT > /dev/null

# Create FCCs header
matlab -nojvm -r "model_file='$MODEL_FILE';" < $HEADER_SCRIPT > /dev/null

HEADER=$(echo -en "Conc_set\tStable_set\tReaction\t"; \
cat /tmp/cbb_reaction_header.long.txt | grep -v "^BioMass_" | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Concatenate FCC data
(echo -e "$HEADER"; find /tmp/skm/fcc/ -name *.tab | xargs -n 32 -P 8 cat) > \
$PARA_DIR/concset_stabstate_rxn_FCCs.tab

# Gzip and retire
pigz $PARA_DIR/concset_stabstate_rxn_FCCs.tab
retire $PARA_DIR/concset_stabstate_rxn_FCCs.tab.gz

# Report step done
echo -e "\n\e[92mStep 2: Done.\e[0m\n"

### 3. EXTRACT SAMPLED PARAMETERS ##############################################

# Report progress
echo -e "\n\e[94mStep 3: Extracting parameters...\e[0m\n"

# Create temporary output directory
mkdir -p $TMP_DIR/par

# Extract parameters
m_vars="indir='$PARA_DIR'; outdir='/tmp/skm/par/'"
matlab -nojvm -r "$m_vars" < $PAREX_SCRIPT > /dev/null

# Store parameters header (bulk created by PAREX_SCRIPT above)
HEADER=$(echo -en "Conc_set\tStable\t"; cat /tmp/skm/par_header.long.txt | \
tr "\n" "\t" | tr -d " " | sed -e 's/\t$/\n/')

# Concatenate parameters data
(echo -e "$HEADER"
find /tmp/skm/par/ -name *.tab | while read File; do
  Conc_set=`echo $File | rev | cut -f 1 -d \/ | rev | cut -f 1 -d \.`
  cat $File | sed -e "s/^/${Conc_set}\t/"
done) > $PARA_DIR/concset_stability_parameters.tab

# Gzip and retire
pigz $PARA_DIR/concset_stability_parameters.tab
retire $PARA_DIR/concset_stability_parameters.tab.gz

# Report step done
echo -e "\n\e[92mStep 3: Done.\e[0m\n"

### 4. PLOT METABOLITE CONCENTRATIONS AND RATIOS VS STABILITY ##################

# Plot with R
$CSTAB_SCRIPT $CONC_FILE $PARA_DIR/met_set_vs_percent_steady.tab

### 5. PLOT FCC HEATMAP ########################################################

# Plot with R
$HTMAP_SCRIPT $PARA_DIR/concset_stabstate_rxn_FCCs.tab.gz

### 6. PLOT KM VS CONCENTRATION ################################################

### 7. DIMENSIONALITY REDUCTION OF KM OVER CONCENTRATION #######################

### 8. CLUSTERING OF KM OVER CONCENTRATION #####################################

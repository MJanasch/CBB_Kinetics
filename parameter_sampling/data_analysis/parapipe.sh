#!/usr/bin/env bash

# Pipeline for extracting and plotting SKM parameter sets

source $1

### 1. CREATE LIST WITH PERCENT STABLE STEADY STATES PER CONCENTRATION SET #####
PARA_DIR="/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2017-09-08/"
STAB_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/percent_steady.m"

matlab -nojvm -r "indir='$PARA_DIR'" < $STAB_SCRIPT > /dev/null

### 2. EXTRACT FCC VALUES ######################################################

# Clean up temporary output directory
TMP_DIR="/tmp/skm"
if [ -d $TMP_DIR ]; then
  rm -r $TMP_DIR
fi

# Create temporary output directory
mkdir $TMP_DIR
mkdir ${TMP_DIR}/skm

# Extract FCCs


# Create FCCs header

# Concatenate and gzip FCC data

# Retire file

### 3. EXTRACT SAMPLED PARAMETERS ##############################################

# Clean up temporary output directory

# Extract

# Create parameters header

# Concatenate and gzip parameters data

# Retire file

### 4. PLOT METABOLITE CONCENTRATIONS AND RATIOS VS STABILITY ##################

### 5. PLOT FCC HEATMAP ########################################################

### 6. PLOT KM VS CONCENTRATION ################################################

### 7. DIMENSIONALITY REDUCTION OF KM OVER CONCENTRATION #######################

### 8. CLUSTERING OF KM OVER CONCENTRATION #####################################

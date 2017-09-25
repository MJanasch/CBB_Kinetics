# Pipeline for extracting and plotting SKM metabolite concentrations
# Input variables

### INPUT FILES ################################################################
MODEL_FILE="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/Data/CBB_Data_170915.mat"
CONC_DIR="/ssd/common/proj/Kinetic_Model/Metabolite_Sampling/Results/2017-09-18"

### PATHS TO SCRIPTS ###########################################################
HEADER_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/metabolite_sampling/data_analysis/extract_header.m"
CONC_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/metabolite_sampling/data_analysis/extract_concentrations.m"
PLOT_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/metabolite_sampling/data_analysis/examine_sampled_concs.R"

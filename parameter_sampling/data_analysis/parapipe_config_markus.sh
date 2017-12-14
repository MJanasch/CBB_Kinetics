# Pipeline for extracting and plotting SKM parameter sets
# Input variables

### INPUT FILES/DIRECTORIES ####################################################
PARA_DIR="/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2017-09-08/"
MODEL_FILE="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/Data/CBB_Data_170908.mat"
CONC_FILE="/ssd/common/proj/Kinetic_Model/Metabolite_Sampling/Results/2017-09-18/all_metabolite_concentrations.tab"

### TEMPORARY DIRECTORY ########################################################
TMP_DIR="/tmp/skm"

### PATHS TO SCRIPTS ###########################################################
STAB_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/percent_steady.m"
FCCEX_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_FCCs.m"
HEADER_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_reaction_header.m"
PAREX_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_parameters.m"
CSTAB_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/concs_vs_stability.R"
HTMAP_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/plot_FCCs_heatmap.R"
MET_NAMES_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/data_analysis/extract_metabolite_header.m"
CCCEX_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/data_analysis/extract_CCCs.m"
HTMAP_CCC_SCRIPT="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/data_analysis/plot_CCCs_heatmap.R"
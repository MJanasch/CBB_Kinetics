# Pipeline for extracting and plotting SKM parameter sets
# Input variables

### INPUT FILES/DIRECTORIES ####################################################
PARA_DIR="/ssd/common/proj/Kinetic_Model/SKM_Sampling/Results/2018-03-10/a/"
MODEL_FILE="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/Data/CBB_Data_180215.mat"
CONC_FILE="/ssd/common/proj/Kinetic_Model/Metabolite_Sampling/Results/2018-03-07/a/all_metabolite_concentrations.tab"
META_FILE="/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/data_analysis/concentrations_and_cofactor_ratios_from_2016.tab"

### TEMPORARY DIRECTORY ########################################################
TMP_DIR="/tmp/skm"

### PATHS TO SCRIPTS ###########################################################
STAB_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/percent_steady.m"
FCCEX_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_FCCs.m"
HEADER_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_reaction_header.m"
PAREX_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/extract_parameters.m"
CSTAB_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/concs_vs_stability.R"
HTMAP_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/plot_FCCs_heatmap.R"
CLSTR_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/cluster_FCCs.R"
CONKM_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/plot_conc_over_Km.R"
CONKAI_SCRIPT="/home/johannes/proj/kimo/tools/maks/parameter_sampling/data_analysis/plot_conc_over_Kai.R"

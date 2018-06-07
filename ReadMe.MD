# Metabolite and Parameter Sampling Scripts for Janasch et al., 2018 (Journal of Experimental Botany)
"Kinetic modeling of the Calvin cycle identifies flux control and stable metabolomes in Synechocystis carbon fixation"

---

### 1. Metabolite Sampling

##### Input data used:

`metabolite_sampling/CBB_NET_Results_XFPK_171026.txt` as `InputNET`- text file with metabolite concentration ranges in [mM]
`parameter_sampling/Data/CBB_Data_XFPK_180307.mat` as `InputDataStructure` - .mat-file containing model structure, equilibrium constants for relevant reactions, direction as they are employed in the model
`NrSampling` - number of sampling rounds, 3.125E6 x 16 employed


##### Metabolite sampling scripts:

`metabolite_sampling/MJanasch_CBB_Metabolite_Sampling.m` - function to sample metabolite concentrations and evaluate thermodynamic feasibility
`metabolite_sampling/MJanasch_CBB_ScriptToCallSamplingFunction.m` - script for calling metabolite sampling function (16x)

##### Output:

16x `MetConcDataSet`, containing all feasible metbolite concentration sets from each of the 16 samplings (3.125E6 each)

---

### 2. Data Analysis of Metabolite Sampling

`metabolite_sampling/data_analysis/metapipe_config` - define directory that contain metabolite sampling output and which model file was used
`metabolite_sampling/data_analysis/metapipe` - pipeline for metabolite sampling analysis

##### Output:

`all_metabolite_concentrations.all_vs_all.png` - metabolite concentration correlations
`all_metabolite_concentrations.tab` - table containing all metabolite concentration sets in mM

---

### 3. Parameter Sampling

##### Input data used:

`.../all_metabolite_concentrations.tab` as `MetConcSamplingData`- .tab-file from metabolite sampling and subsequent data analysis
`parameter_sampling/Data/CBB_Data_XFPK_180307.mat` as `InputDataStructure` - .mat-file containing model structure, equilibrium constants for relevant reactions, direction as they are employed in the model
`iterations` - number of parameter sampling rounds per metabolite concentration set, here: 1000
`NrOfMetDataSet` - number of concentration set 

##### Parameter sampling scripts:

`parameter_sampling/MJanasch_SKM_Sampling_Code.m` - function to sample parameter values around metabolite concentrations, number of samplings defined as `ìterations`
`parameter_sampling/MJanasch_CBB_ScriptToCallSKMSamplingFunction.m` - script for calling parameter sampling function for each concentration set 

##### Output:

 `DataOut` for each metabolite concentration set, containing the concentrations used, 1000 (`iterations`) parameter sets, their stability, and Flux Control Coefficients ('CJ_rec'), amongst other information

---

### 4. Data Analysis of Parameter Sampling

`parameter_sampling/data_analysis/parapipe_config` - define directory that contain metabolite sampling output, parameter sampling output and which model file was used
`parameter_sampling/data_analysis/parapipe` - pipeline for parameter sampling analysis

##### Output:
(amongst others)

`met_set_vs_percent_steady.tab` - table containing percentage of parameter sets resulting in stable steady states for each concentration set
`concs_vs_stability.pdf` - tendencies for all metabolites causing most (top 10%) and fewest (bottom 10%) stable steady states
`conc_over_Km.pdf` - tendencies for saturation states for each metabolite-reaction-pair to cause stable or unstable steady states
`FCCs_heatmap.pdf` - Heatmap of the Flux Control Coefficients (FCCs) for each reaction pair
`effector_patterns_pvclust.pdf` - Clustering of FCCs according to their ability of controlling other reactions
`target_patterns_pvclust.pdf` - Clustering of FCCs according to their ability of being controlled

---

## Dependencies
- MatLab ≥ v. 2015a
- SBML Toolbox for MatLab ≥ 4.1.0 ([http://sbml.org/Software/SBMLToolbox](http://sbml.org/Software/SBMLToolbox))

- R ≥ 3.0 with following libraries
	data.table
	dplyr
	GGally
	ggplot2
	ggrepel
	pvclust
	reshape2
	scales

---

## Author
Markus Janasch (<markus.janasch@scilifelab.se>)
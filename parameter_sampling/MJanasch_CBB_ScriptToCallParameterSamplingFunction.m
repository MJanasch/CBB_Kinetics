% This script calls the fast parameter sampling and coefficient calculating 
% function for the CBB Model, adapted from Murabito et al., 2014 (Steuer)
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-10-12, last modified: -

%% Add SBML-folders to path on Monty
addpath('/hdd/common/tools/sbml/libSBML-5.15.0-matlab/')
addpath('/hdd/common/tools/sbml/SBMLToolbox-4.1.0/')

%% Add path to the SKM/Parameter Sampling function

addpath('/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/');

MetConcData_RAW = importdata(MetConcSamplingData);



%% Call SKM-algorithm
[DataOut] = MJanasch_Parameter_Sampling_Code(Iterations,InputDataStructure,MetConcData_RAW.data(NrOfMetDataSet,:),MetConcData_RAW.textdata);


save(DataSetOut,'DataOut');
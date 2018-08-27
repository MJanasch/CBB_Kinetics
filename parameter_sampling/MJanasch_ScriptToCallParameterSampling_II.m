% This script calls the parameter sampling and coefficient calculating 
% functions for the kinetic models
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-03-22, last modified: 2018-08-26

%% Add path to the Parameter Sampling function

addpath('/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/');

MetConcData_RAW = importdata(MetConcSamplingData);

%MJanasch_CreateFunction_CalDFODC(InputDataStructure)

[DataOut] = MJanasch_Parameter_Sampling_II(Iterations,InputDataStructure,MetConcData_RAW.data(NrOfMetDataSet,:),MetConcData_RAW.textdata);
save(DataSetOut,'DataOut');
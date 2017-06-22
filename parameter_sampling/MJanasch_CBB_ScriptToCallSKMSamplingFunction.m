% This script calls the SKM parameter sampling and coefficient calculating 
% function of Murabito et al., 2014 (Steuer)
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-11, last modified: 2017-05-12

%% Add SBML-folders to path on Monty
addpath('/hdd/common/tools/sbml/libSBML-5.15.0-matlab/')
addpath('/hdd/common/tools/sbml/SBMLToolbox-4.1.0/')

MetConcData_RAW = importdata('Data/metabolite_concs_16x5e8_samples_2017-05-10.tab');




%% Call SKM-algorithm
[DataOut] = MJanasch_SKM_Sampling_Code(Iterations,MetConcData_RAW.data(NrOfMetDataSet,:));


save(DataSetOut,'DataOut');
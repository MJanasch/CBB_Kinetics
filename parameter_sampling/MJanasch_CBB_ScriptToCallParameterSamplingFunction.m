% This script calls the fast parameter sampling and coefficient calculating 
% function for the CBB Model, adapted from Murabito et al., 2014 (Steuer)
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-10-12, last modified: -

%% Add SBML-folders to path on Monty
addpath('/hdd/common/tools/sbml/libSBML-5.15.0-matlab/')
addpath('/hdd/common/tools/sbml/SBMLToolbox-4.1.0/')

%% Add path to the Parameter Sampling function

addpath('/ssd/common/proj/Kinetic_Model/maks/parameter_sampling/');

MetConcData_RAW = importdata(MetConcSamplingData);

if exist('ModelType')
    if ~strncmp(ModelType,'CBB',3) && ~strncmp(ModelType,'XFPK',4)
        fprintf('%s\n', 'ModelType not specified. Choose either CBB or XFPK');
    else  
        [DataOut] = MJanasch_Parameter_Sampling_Code(Iterations,InputDataStructure,MetConcData_RAW.data(NrOfMetDataSet,:),MetConcData_RAW.textdata,ModelType);
        save(DataSetOut,'DataOut');
    end
else
    fprintf('%s\n', 'ModelType not specified. Choose either CBB or XFPK');
end

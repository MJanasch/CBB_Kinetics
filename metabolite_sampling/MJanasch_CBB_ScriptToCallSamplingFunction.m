% This script calls the Metabolite-Sampling function
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-08, last modified: 2017-07-04

%% Add path to the Metabolite Sampling function

addpath('/ssd/common/proj/Kinetic_Model/maks/metabolite_sampling/');

%% Actual Script

% Seed = 0;     % Default value for Seed
rng(Seed)       % Define Seed for random number generator "rand", has to be
                % set outside matlab for running on several cores
%[Y,MetConcDataSet,Infeasible_Reactions] = MJanasch_CBB_Metabolite_Sampling(NrSampling,InputDataStructure,InputNET);
[Y,MetConcDataSet,dGDataSet] = MJanasch_CBB_Metabolite_Sampling(NrSampling,InputDataStructure,InputNET);


save(outfile,'MetConcDataSet');
%save(outfile_2,'Infeasible_Reactions');

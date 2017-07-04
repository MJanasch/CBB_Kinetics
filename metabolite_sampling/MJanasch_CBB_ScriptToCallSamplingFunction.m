% This script calls the Metabolite-Sampling function
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-08, last modified: 2017-07-04


% Seed = 0;     % Default value for Seed
rng(Seed)       % Define Seed for random number generator "rand", has to be
                % set outside matlab for running on several cores
[Y,MetConcDataSet] = MJanasch_CBB_Metabolite_Sampling(NrSampling,InputDataStructure,InputNET);


save(outfile,'MetConcDataSet');
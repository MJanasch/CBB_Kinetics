% This script calls the Metabolite-Samnpling function
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-05-08, last modified: 2017-05-10


rng(Seed)
[Y,MetConcDataSet] = MJanasch_CBB_Metabolite_Sampling(NrSampling);


save(outfile,'MetConcDataSet');
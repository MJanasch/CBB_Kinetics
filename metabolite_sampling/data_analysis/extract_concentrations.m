% Extract concentrations from a set of mat files

% Create list of infiles; assuming only one set of sampled concentrations
% is present in the input directory
infiles = dir(fullfile(indir,'*.mat'));
infiles = {infiles.name}.';

% Iterate over the infiles
n = 1;
load(char(fullfile(indir, infiles(n))));
concentrations = MetConcDataSet;

for n = 2:length(infiles)
  % Load data
  load(char(fullfile(indir, infiles(n))));
  % Add concentration data to end of concentrations matrix
  concentrations = horzcat(concentrations, MetConcDataSet);
end

% Write concentration matrix to file
dlmwrite('/tmp/cbb_metabolite_concs.tab', transpose(concentrations), 'delimiter', '\t', 'precision', '%10d');

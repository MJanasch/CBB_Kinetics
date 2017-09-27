% Extract percent stable steady states from a set of CBB SKM mat files

% Create list of infiles; assuming only one set of sampled parameters
% is present in the input directory
infiles = dir(fullfile(indir,'*.mat'));
infiles = {infiles.name}.';
N = length(infiles);

% Iterate over the infiles
p_stab = [];

for n = 1:N
  fprintf(2, '%3.1f%%\r', n/N*100)
  % Load data
  load(char(fullfile(indir, infiles(n))));
  % Acquire metabolite concentration set ID
  infile_split = strsplit(char(infiles(n)), '_');
  infile_n = str2num(char(infile_split(3)));
  % Add stability data to end of stability matrix
  n_stable = sum(DataOut.StabilityIndicator);
  n_tot = length(DataOut.StabilityIndicator);
  p_stab = [p_stab; [infile_n, n_stable/n_tot*100]];
end

fprintf(2, '%3.1f%%\n', n/N*100)

% Write concentration matrix to file
outfile = 'met_set_vs_percent_steady.tab';
dlmwrite(fullfile(indir, outfile), p_stab, 'delimiter', '\t');

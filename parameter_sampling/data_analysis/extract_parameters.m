% Extract parameters from a set of CBB SKM mat files

% Create list of infiles; assuming only one set of sampled parameters
% is present in the input directory
infiles = dir(fullfile(indir,'*.mat'));
infiles = {infiles.name}.';
N = length(infiles);

% Iterate over the infiles
for n = 1:N
  fprintf(2, '%3.1f%%\r', n/N*100)
  % Load data
  load(char(fullfile(indir, infiles(n))));
  % Acquire metabolite concentration set ID
  infile_split = strsplit(char(infiles(n)), '_');
  infile_n = str2num(char(infile_split(3)));
  % Extract parameters and vector indicating stability or not
  parameters = horzcat(DataOut.StabilityIndicator, DataOut.Parameters);
  % Write to tab-delimited file
  outfile = fullfile(outdir, strcat(num2str(infile_n), '.tab'));
  dlmwrite(outfile, parameters,'delimiter', '\t');
end

% Save a header
dlmwrite('/tmp/skm/par_header.long.txt', char(DataOut.ParID), 'delimiter', '')

fprintf(2, '%3.1f%%\n', n/N*100)

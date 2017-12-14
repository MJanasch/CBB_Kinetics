% Extract concentration control coefficients from a set of CBB SKM mat files

% Prepare column with metabolite numbers
load(model_file);
[n_met,n_rxn] = size(SFull);
met_col = transpose(1:n_met);

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
  % Prepare column with metabolite sample set numbers
  set_col = ones(n_met,1) * infile_n;
  % Determine size of FCCs matrix
  s = size(DataOut.CS_rec);
  % Determine
  if size(s) == [1 2]
    I = 1;                  % Only one stable steady state
  else
    % Skip if there are no stable steady states
    if s(3) == 0
      continue
    else
      I = s(3);             % Number of stable steady states
    end
  end
  % Iterate over the stable steady states
  for i = 1:I
    i_col = ones(n_met,1) * i;
    % Extract and modify each matrix
    C = DataOut.CS_rec(:,:,i);
    C = [set_col i_col met_col C];
    % Write to tab-delimited file
    outfile = fullfile(outdir, strcat(num2str(infile_n), '_', num2str(i), '.tab'));
    dlmwrite(outfile, C,'delimiter', '\t');
  end
end
fprintf(2, '%3.1f%%\n', n/N*100)

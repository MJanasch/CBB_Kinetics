% Extract the metbolites header from a CBB SKM model file
load(model_file);
Floating_Metabolites = N.species(1:size(SFull,1));
header = char({Floating_Metabolites.name}.');
outfile = fullfile(outdir, 'cbb_met_list.long.txt');
dlmwrite(outfile, header, 'delimiter', '');

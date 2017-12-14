% Extract the metbolites header from a CBB SKM model file
load(model_file);
header = char({N.species.name}.');
outfile = fullfile(outdir, 'cbb_met_list.long.txt');
dlmwrite(outfile, header, 'delimiter', '');

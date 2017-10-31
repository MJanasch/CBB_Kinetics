% Extract the reactions header from a CBB SKM model file
load(model_file);
header = char({N.reaction.id}.');
outfile = fullfile(outdir, 'cbb_reaction_header.long.txt');
dlmwrite(outfile, header, 'delimiter', '');

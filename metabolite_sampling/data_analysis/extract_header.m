% Extract the metabolites header from a CBB SKM model file
load(infile);
header = char({N.species.id}.');
dlmwrite('/tmp/cbb_metabolite_header.long.txt', header, 'delimiter', '');

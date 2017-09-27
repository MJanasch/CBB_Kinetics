% Extract the reactions header from a CBB SKM model file
load(model_file);
header = char({N.reaction.id}.');
dlmwrite('/tmp/cbb_reaction_header.long.txt', header, 'delimiter', '');

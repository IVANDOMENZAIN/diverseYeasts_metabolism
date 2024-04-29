function mutant = getGALmutant(model)
galGenes = {'Seq_1935'	'Seq_4294'	'Seq_3460'	'Seq_2479'	'Seq_3332'};
mutant   = removeGenes(model,galGenes,false,false,true);
end
%strategy is to first open the model, set the genes in a file (variable)
%create a file (variable) for the orthologues
%compare the range from the file which contains the 131 genes and run them
%against the model.genes file. The ones that do not match, to be appended
%to the end of the model.genes? dont know if that is even possible


%load the model
load('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/models/candida_intermedia/cintGEM_oxido.mat');
model.genes = strrep(model.genes,'Candida_intermedia@','');
%import the orthologue_unassigned genes file
Y = readtable('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups_UnassignedGenes.txt','delimiter','\t');
%Checking the unassigned genes vs the genes in the model
unassigned_genes = Y.model_Cint(~cellfun(@isempty,Y.model_Cint));
C = find(ismember(model.genes,unassigned_genes));

%checking the presence of the unassigned genes in the orthogroups
orthogroups = readtable('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups.txt','delimiter','\t');
[presence,idxs] = ismember(model.genes,orthogroups.model_Cint);
% it works!
idxs2 = find(presence);
idxs = idxs(presence);
model.orthologues = model.genes;
model.orthologues(idxs2) = orthogroups.data_Cint(idxs);

save('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/models/candida_intermedia/cintGEM_oxido.mat','model');
%strategy is to first open the model, set the genes in a file (variable)
%create a file (variable) for the orthologues
%compare the range from the file which contains the 131 genes and run them
%against the model.genes file. The ones that do not match, to be appended
%to the end of the model.genes? dont know if that is even possible


%load the model
load('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/models/candida_intermedia/cintGEM_oxido.mat');
%import the orthologue_unassigned genes file
Y = tdfread('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups_UnassignedGenes.tsv','tab');
% I just mix both
C = union(model.genes, Y.model_Cint);
%will just check if anything is repeating
model.genes = unique(C);


%checking how unique actually works
load('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/models/candida_intermedia/cintGEM_oxido.mat');
%import the orthologue_unassigned genes file

Y = tdfread('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups_UnassignedGenes.tsv','tab');

T = Y.model_Cint;
dude_before = []
for i = 1:length(T)
    dude_before = strcmp(model.genes,T(i));
end
        
C = union(model.genes, Y.model_Cint);
dude_after = []
for i = 1:length(T)
    dude_after = strcmp(C,T(i));
end
        
%repMets_analysis

%load model 
load('../models/candida_intermedia/cintGEM_oxido.mat')
cSources = {'gal' 'lac' 'cel' 'xyl'};
mkdir('../results/reporter_metabolites')
for cSource=cSources
    %load data
    str = cSource{1};
    disp(str)
    disp(' ')
    DEresults = readtable(['../results/RNA_DE_analysis/RNA_DE_glu_vs_' str '.txt'],'delimiter','\t');
    DEmapping = readtable(['../results/RNA_DE_analysis/RNA_2_model_glu_vs_' str '.txt'],'delimiter','\t');
    
    %prepare inputs for the rep mets function
    genes           = DEmapping.modelGenes;
    genePvalues     = DEmapping.adjPval;
    geneFoldChanges = DEmapping.log2FC;
    
    %Run reporter metabolites analysis!
    outputFile      = ['../results/reporter_metabolites/repMets_glu_vs_' str '.txt'];
    repMets=reporterMetabolites(model,genes,genePvalues,true,outputFile,geneFoldChanges);
end
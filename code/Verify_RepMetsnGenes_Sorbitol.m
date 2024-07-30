%explore_CintGEM_wFSEOF
load('../models/candida_intermedia/cintGEM_oxido.mat')
%identify relevant rxns and mets associated to L-glucitol in the model 
posGluc = find(contains(model.metNames,'L-glucitol'));

%Verify glucitol presence 
disp(model.metNames(posGluc))
disp(model.metComps(posGluc))
disp(model.compNames)

%the model shows presence of glucitol in the cytoplasm and excretes it
%outisde as well. This does not confirm glucitol associated to the
%oxidopathway yet
posGluc = posGluc(1);
GlucRxns = find(model.S(posGluc,:));

%now we found two reactions associted with the cytoplasmic sorbitol. Time
%to check the reactions.
constructEquations(model,GlucRxns,true)

%the model has both conversion of sorbitol to sorbose but also
%extracellular transport of sorbitol. The conversion from sorbitol to
%sorbose is a step in the oxido pathway. 

GlucRxns =GlucRxns(1);
GlucGenes =model.grRules(GlucRxns);
disp(GlucGenes)

%RepMet results showed three connected genes Seq_1183,Seq_2189,Seq_4936 but
%here we see only one. Lets explore the other genes

otherGenes = {'Seq_2189', 'Seq_1183', 'Seq_4936'};

%Here we go gene by gene finding associated reacitons and print the
%formulas
for i=1:length(otherGenes)
    genePos = find(strcmp(model.genes,otherGenes{i}));
    geneAssRxns = find(model.rxnGeneMat(:,genePos));
    disp(otherGenes{i})
    constructEquations(model,geneAssRxns,true)
    disp(' ')
end

%we found that Seq_2189 is annotated as XYL2 in the genome and in the RNA_2_model_glu_vs_gal file is 5 times log2 fold upregulated. the p value is significant.  





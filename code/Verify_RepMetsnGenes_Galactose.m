%We found that in the reporter metabolite analysis with glucose vs lactose,
%galactose was the rep metabolite. so, similar to whats done with sorbitol,
%I am checking the genes associated with galactose in the model.
load('../models/candida_intermedia/cintGEM_oxido.mat')
%identify relevant rxns and mets associated to D-galactose in the model 
posGal = find(strcmp(model.metNames,'D-galactose'));
%Verify galactose presence 
disp(model.metNames(posGal))
disp(model.metComps(posGal))
disp(model.compNames)

%the model shows presence of galactose in the cytoplasm, vacuole and excretes it
%outisde as well. This does not confirm galactose associated to the
%oxidopathway yet
%then lets check for the 1st D-galactose, the rxns associated to it. 
posGal = posGal(1);
GalRxns = find(model.S(posGal,:));

%now we found 8 reactions associted with the cytoplasmic galactose. Time
%to check the reactions.
constructEquations(model,GalRxns,true)

% we see mainly D-galactose in association with leloir and oxido
% pathway.This means that D-galactose to galactitol is also a
% transcriptionally important reaction.

% looking for the genes in the model that are associated to these reactions
newGalRxns = GalRxns(8);
GalGenes = model.grRules(newGalRxns);
disp(GalGenes)
%These are genes YBR020W(Seq_1935) which is not relevant based on RNAseq expression
%data. And...YBR020W(Seq_4294)

%similarly for the second instance of D-galactose that is associted with
%transport reactions
posGal2 = posGal(2);
GalRxns2 = find(model.S(posGal2,:));
constructEquations(model,GalRxns2,true)
GalRxns2 = GalRxns2(1);
GalGenes2 = model.grRules(GalRxns2);
disp(GalGenes2)
%the only gene that is relevant is STL1(Seq_847) which is downregulated and
%has high p value. No idea what this is. But its associated to D-galactose
%transport, so not sure. 
%lets check the third instance of galactose

%RepMet results showed three connected genes Seq_1183,Seq_2189,Seq_4936 but
%here we see only one. Lets explore the other genes
posGal3 = posGal(3);
GalRxns3 = find(model.S(posGal3,:));
constructEquations(model,GalRxns3,true)
GalRxns3 = GalRxns3(1);
GalGenes3 = model.grRules(GalRxns3);
disp(GalGenes3)


otherGenes = {'Seq_4642','Seq_1183','Seq_2189','Seq_1160','Seq_4006','Seq_758','Seq_5357','Seq_2923','Seq_5493','Seq_4936','Seq_4728'};

%Here we go gene by gene finding associated reacitons and print the
%formulas
for i=1:length(otherGenes)
    genePos = find(strcmp(model.genes,otherGenes{i}));
    geneAssRxns = find(model.rxnGeneMat(:,genePos));
    disp(otherGenes{i})
    constructEquations(model,geneAssRxns,true)
    disp(' ')
end

%we found that Seq_2189 is annotated as XYL2 in the genome and in the RNA_2_model_glu_vs_gal file is 5 times log2 fold upregulated and 2.5 times in lactose. the p value is significant.
%found also that Glucitol is present once again in the list of reporter
%metabolites. Checked the associated genes (Seq_5357 = YALI0B16192G is 3
%fold upregulation with low p-value), (Seq_2923 = SOU2 -1.5 fold change and
%low p-value), (Seq_2552 = lad highly upreg in galactose)
%also found galactose in the repmets_glu_gal but no surprises. another one
%is Seq_2272 = YMR226c which I have deleted, but its not really
%significant in lactose, but is in galactose.






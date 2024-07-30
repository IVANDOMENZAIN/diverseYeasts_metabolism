%An orthofinder run comparing a T. reesei LXR4 sequence vs. the FASTA file
%used for the model reconstruction revealed that this gene has an
%orthologue in C. intermedia (Seq_2272)
clc
clear
%Let's try to find the gene in our model
load('../../models/candida_intermedia/cintGEM_oxido_orthologs.mat')
%search gene in model
genePos = find(strcmpi(model.genes,'Seq_2272'));
%The gene is present let's see if it has been detected also in our RNA
%dataset
ortholog = model.orthologues(genePos);
% YES!!! the gene has also an associated ID in the RNAseq datasets
%let's find this gene in the RNAseq dataset, what's its current annotation?
dataset  = readtable('../../data/RNAseq/normalized_counts.txt','delimiter','\t');
lxr4pos  = find(strcmpi(dataset.genes,ortholog));
lxr4data = dataset(lxr4pos,:);
%Next step: search for our gene in DE results
conditions = {'glu' 'gal' 'lac'};
newTable = table();
for j = 1:length(conditions)
    cond = conditions{j};
    if ~strcmpi(cond,'glu')
        DE_results = readtable(['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' cond '.txt'],'delimiter','\t');
        index = find(strcmpi(DE_results.Row,ortholog));
        newTable = [newTable; [DE_results(index,:), conditions(j)]];
    end
end
%It was found that our gene is significantly upreg in presence of galactose
%let's search for the rest of the pathway
genes = {'xyl1' 'xyl1_2' 'GAL1' 'lad' 'XYL2'};
for i = 1:length(genes)
    gene = genes{i};
    idx  = find(strcmpi(dataset.geneNames,gene));
    geneDataID = dataset.genes(idx);
    %Apend the DE results for these genes to our new table
    for j = 1:length(conditions)
        cond = conditions{j};
        if ~strcmpi(cond,'glu')
            DE_results = readtable(['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' cond '.txt'],'delimiter','\t');
            index = find(strcmpi(DE_results.Row,geneDataID));
            newTable = [newTable; [DE_results(index,:), conditions(j)]];
        end
    end
    
end
newTable = sortrows(newTable,'adjPVal','ascend');
newTable = sortrows(newTable,'Var8','ascend');
%Save results
writetable(newTable,'../../results/RNA_DE_analysis/DE_oxido_leloir_genes_lac_gal.txt','QuoteStrings',false,'Delimiter','\t')

%going back to the model: what's up with lxr4?
lxr4Rxns = find(contains(model.grRules,'Seq_2272'));
%The genes exists in the model and encodes for what has been reported on
%uniprot for its cerevisiae orthologues.
%Let's correct gene association for the lxr4 rxns (oxi/red pathway) in the model
try
    idx = find(contains(model.orthologues,'lxr4'));
    idx = find(contains(model.grRules,'lxr4'));
    %correct gene association 
    model.grRules{idx} = 'Candida_intermedia@Seq_2272';
    model.rxnGeneMat(idx,:) = 0*model.rxnGeneMat(idx,:);
    model.rxnGeneMat(idx,genePos) = 1;
catch
    disp('G0RNA2 was not found in the model genes')
end
[grRules, rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
model.genes = model.genes(1:(end-1));
%Update model 
model.orthologues = model.orthologues(1:(end-1));
model.proteins = model.proteins(1:(end-1));
idx = find(strcmpi(model.genes,'Candida_intermedia@Seq_2272'));
model.geneShortNames(idx) = {'LXR4'};
model.proteins(idx) = {'lxr4'};
model.geneShortNames = model.geneShortNames(1:(end-1));
%Correct grRules and rxnGEneMat
model.grRules = strrep(model.grRules,'Candida_intermedia@','');
[grRules,rxnGeneMat] = standardizeGrRules(model);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;

%generate version-controllable files
formulas = constructEquations(model);
rxns = model.rxns;
rxnNames = model.rxnNames;
grRules = model.grRules;
modelTable = table(rxns,rxnNames,formulas, grRules);
writetable(modelTable,'../../models/candida_intermedia/cintGEM_oxido_orthologs_curated.txt','WriteVariableNames',true,'Delimiter','\t','QuoteStrings',false);

%add version control
genes = model.genes;
model.geneShortNames = strrep(model.geneShortNames,'Candida_intermedia@','');
shortnames = model.geneShortNames;
orthologues = model.orthologues;
proteins = model.proteins;
gene_table = table(genes,shortnames,orthologues,proteins);
writetable(gene_table,'../../models/candida_intermedia/gene_table_CintOxido_orthologues_curated.txt','Delimiter','\t','QuoteStrings',false);

save('../../models/candida_intermedia/cintGEM_oxido_orthologs_curated.mat','model')



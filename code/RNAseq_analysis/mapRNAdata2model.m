%Load model 
load('../../models/candida_intermedia/cintGEM.mat')
%Correct grRules field
%model.grRules = strrep(model.grRules,'Candida_intermedia@','');
%[grRules,rxnGeneMat] = standardizeGrRules(model);
%model.grRules = grRules;
%model.rxnGeneMat = rxnGeneMat;
threshold = 0.01;
%Iterate through each condition
conditions = {'gal' 'lac' 'cel' 'xyl'};
nanMat   = zeros(length(model.orthologues),5);
nanMat(:,:) = nan;

orthoIDs = model.orthologues;
for k=1:length(orthoIDs)
    if contains(model.orthologues{k},';')
        x = strsplit(model.orthologues{k},';');
        orthoIDs{k} = x{1};
    end 
end

for j = 1:length(conditions)
    newTable = table();
    cond = conditions{j};
    disp(cond)
    %Get RNAseq like gene IDs from model
    newTable.dataGenes = (model.orthologues);
    %open RNAseq results file for the given condition 
    DE_results = readtable(['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' cond '.txt'],'delimiter','\t');
    %find which genes in the model are also part of the RNAseq dataset


    [iA,iB] = ismember(orthoIDs,DE_results.Row);
    newTable.geneNames = cell(height(newTable),1);
    newTable.geneNames(find(iA)) = DE_results.geneNames(iB(iB>0));
    newTable.modelGenes = model.genes;
    newTable.proteins = model.proteins;
    
    newTable.counts_ref = zeros(height(newTable),1);
    newTable.counts_ref(find(iA)) = DE_results.Mean_Ref(iB(iB>0));
    
    newTable.counts_Csrc = zeros(height(newTable),1);
    newTable.counts_Csrc(find(iA)) = DE_results.Mean__Csrc(iB(iB>0));    
    
    newTable.log2FC = zeros(height(newTable),1);
    newTable.log2FC(find(iA)) = DE_results.Log2_FC(iB(iB>0));    

    newTable.adjPval = zeros(height(newTable),1);
    newTable.adjPval(find(iA)) = DE_results.adjPVal(iB(iB>0));  
    %Map the rxns and mets that are linked to each gene
    %Create a metGeneMatrix
    metGeneMat = logical(model.S)*(model.rxnGeneMat);
    newTable.rxns = cell(height(newTable),1);
    newTable.mets = cell(height(newTable),1);
    for i=1:length(model.genes)
       rxns = find(model.rxnGeneMat(:,i));
       mets = find(metGeneMat(:,i));
       newTable.rxns{i} =rxns;
       newTable.mets{i} =mets;
    end
    DEtable = newTable(newTable.adjPval <= threshold & newTable.counts_ref>0 & newTable.counts_Csrc>0,:);
    upReg = sum(DEtable.log2FC>0);
    dReg  = sum(DEtable.log2FC<0);
    disp(['There are ' num2str(height(DEtable)) ' DE genes, from which ' num2str(upReg) ' are upregulated, and ' num2str(dReg) ' are downregulated'])
	writetable(newTable,['../../results/RNA_DE_analysis/RNA_2_model_glu_vs_' cond '.txt'],'delimiter','\t','QuoteStrings',false,'WriteRowNames',true)

end
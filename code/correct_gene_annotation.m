% orthogroups     = readtable('../orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups.txt','delimiter','\t');
% newDataset = getFastaIDs;

load('../models/candida_intermedia/cintGEM_curated.mat')
model.orthologues(1071) = {'galactitol_dh'};
t = table(model.genes,model.geneShortNames,model.orthologues);
writetable(t,'../results/modelGenes.txt','delimiter','\t','QuoteStrings',false)
formulas = constructEquations(model);
t = table(model.rxns,model.rxnNames,formulas,model.grRules);
writetable(t,'../results/model_grRules.txt','delimiter','\t','QuoteStrings',false)

inconsistencies = find(contains(model.orthologues,'Seq'));
%there were 154 inconsitencies found 14.39%
anomalies = [];
genes2add = [];
for i=1:numel(model.orthologues)
    seqGene = model.orthologues(i);
    if contains(seqGene,'Seq_')
        position = find(contains(orthogroups.model_Cint,seqGene));
        if ~isempty(position)
            if numel(position)==1
                if ~isempty(orthogroups.data_Cint{position})
                    genesModel = strsplit(orthogroups.model_Cint{position},',');
                    genesRNASQ = strsplit(orthogroups.data_Cint{position},',');
                    genesModel = strtrim(genesModel);
                    genesRNASQ = strtrim(genesRNASQ);
                    if length(genesModel) <= length(genesRNASQ)    
                        for j=1:length(genesModel)
                            b = find(strcmp(model.orthologues,genesModel{j}));
                            if ~isempty(b)
                                x = find(contains(newDataset.IDs_1,genesRNASQ{j}));
                                model.orthologues(b) = newDataset.IDs_2(x);
                            end
                        end
                        a = ~ismember(genesModel,model.orthologues);
                        genes2add  = [genes2add; genesModel(a)'];
                    else
                        for j=1:length(genesRNASQ)
                            x = find(strcmp(model.genes,genesModel(j)));
                            y = find(contains(newDataset.IDs_1,genesRNASQ{j}));
                            model.orthologues(x) = newDataset.IDs_2(y);
                        end
                        %newGenes = length(genesRNASQ) - length(genesModel)
                    end
                end
            else
                anomalies = [anomalies;seqGene];
            end
        else
            %disp(seqGene)
            %disp(orthogroups.model_Cint(position))
        end
    end
end

t = table(model.genes,model.geneShortNames,model.orthologues);
writetable(t,'../results/modelGenes.txt','delimiter','\t','QuoteStrings',false)
formulas = constructEquations(model);
t = table(model.rxns,model.rxnNames,formulas,model.grRules);
writetable(t,'../results/model_grRules.txt','delimiter','\t','QuoteStrings',false)
save('../models/candida_intermedia/cintGEM_curated2.mat','model')

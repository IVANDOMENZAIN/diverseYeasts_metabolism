orthogroups = readtable('../orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups.txt','delimiter','\t');
newDataset  = getFastaIDs;
load('../models/candida_intermedia/cintGEM_oxido_orthologs_curated.mat')
inconsistencies = find(contains(model.orthologues,'Seq'));
%there were 154 inconsitencies found 14.39% of the model's genes
anomalies = [];
genesIDs = [];
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
                            b = find(strcmp(model.genes,genesModel{j}));
                            if ~isempty(b)
                                x = find(contains(newDataset.IDs_1,genesRNASQ{j}));
                                model.orthologues(b) = newDataset.IDs_2(x);
                                last_b = b;
                            end
                        end
                        if length(genesModel) < length(genesRNASQ) 
                            gene = model.genes(last_b);
                            for k=(j+1):length(genesRNASQ) 
                                nextOrtholog = genesRNASQ{j+1};
                                y = find(contains(newDataset.IDs_1,nextOrtholog));
                                genesIDs  = [genesIDs; newDataset.IDs_2(y)];

                                model.orthologues{last_b} = [model.orthologues{last_b} '; ' newDataset.IDs_2{y}];
                                rxns = find(contains(model.grRules,gene));
                                for l=1:length(rxns)
                                    if ~contains(model.grRules(rxns(l)),' and ')
                                        model.grRules{rxns(l)} = [model.grRules{rxns(l)} ' or ' newDataset.IDs_2{y}];
                                    else
                                        disp(model.grRules(rxns(l)))
                                    end
                                end
                            end
                        end
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
genes2add= struct();
genes2add.genes = genesIDs;
genes2add.geneShortNames = genesIDs;
model = addGenesRaven(model,genes2add);
[a,b] = standardizeGrRules(model,false);
model.grRules = a;
model.rxnGeneMat = b;
model.orthologues = [model.orthologues;  genesIDs];
model.proteins = [model.proteins; genesIDs];
genes = [];
%fix inconsistencies manually, these are
gene = {'Seq_127'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ54830.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_278'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ49614.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_398'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ53714.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_135'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ52714.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_80'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ56702.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_118'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ47884.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_194'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ51356.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_327'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ47430.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_292'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ50208.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_594'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ47281.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_5004'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ57507.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_652'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ56088.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_623'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ55088.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_4649'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ50241.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;

gene = {'Seq_4552'};
x = find(contains(orthogroups.model_Cint,gene));
ortholog = {'SGZ50241.1'};
y = find(contains(newDataset.IDs_1,ortholog));
ortholog = newDataset.IDs_2(y);
x = find(strcmp(model.genes,gene));
model.orthologues(x) = ortholog;


[a,b] = standardizeGrRules(model,false);
model.grRules = a;
model.rxnGeneMat = b;
%generate version-controllable files
formulas = constructEquations(model);
rxns = model.rxns;
rxnNames = model.rxnNames;
grRules = model.grRules;
modelTable = table(rxns,rxnNames,formulas, grRules);
writetable(modelTable,'../models/candida_intermedia/cintGEM_curated.txt','WriteVariableNames',true,'Delimiter','\t','QuoteStrings',false);

%add version control
genes = model.genes;
shortnames = model.geneShortNames;
orthologues = model.orthologues;
proteins = model.proteins;
gene_table = table(genes,shortnames,orthologues,proteins);
writetable(gene_table,'../models/candida_intermedia/gene_table_curated.txt','Delimiter','\t','QuoteStrings',false);

save('../models/candida_intermedia/cintGEM_gene_curated.mat','model')

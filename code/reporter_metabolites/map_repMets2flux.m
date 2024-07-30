%repMets_analysis

%load model 
load('../../models/candida_intermedia/cintGEM_curated2.mat')
%Correct model grRules
model.grRules = strrep(model.grRules,'Candida_intermedia@','');
[grRules,rxnGeneMat] = standardizeGrRules(model,false);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;

cSources = {'gal'};
mkdir('../../results/reporter_metabolites')
multiDim_data = table();

%Exclude highly connected Mets
toExclude = {'carbon dioxide' 'ATP' 'ADP' 'AMP' 'phosphate' 'diphosphate' 'H+' ...
             'CMP' 'GDP' 'GTP' 'H2O' 'CTP' 'coenzyme A'};

for cSource=cSources
    %load data
    str = cSource{1};
    disp(str)
    disp(' ')
    DEresults = readtable(['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' str '.txt'],'delimiter','\t');
    DEmapping = readtable(['../../results/RNA_DE_analysis/RNA_2_model_glu_vs_' str '.txt'],'delimiter','\t');
    nonZeros  = sum(DEmapping.counts_ref+DEmapping.counts_Csrc,2)>0;
    DEmapping = DEmapping(nonZeros,:);
    %
    if isempty(multiDim_data)
        multiDim_data.genes = DEmapping.modelGenes;
    end
    eval(['multiDim_data.pVal_' str '=DEmapping.adjPval;'])
        
    %prepare inputs for the rep mets function
    genes           = DEmapping.modelGenes;
    genePvalues     = DEmapping.adjPval;
    geneFoldChanges = DEmapping.log2FC;
    
    %Run reporter metabolites analysis!
    repMets=get_repMets(model,genes,genePvalues,toExclude,true,outputFile,geneFoldChanges);
    metsUP = repMets(2).mets(find(repMets(2).metPValues<0.001));
    metsDN = repMets(3).mets(find(repMets(3).metPValues<0.001));
    DEgUP = genes(genePvalues<0.001 & geneFoldChanges>0);
    DEgDN = genes(genePvalues<0.001 & geneFoldChanges<0);

    [a,b] = ismember(metsUP,model.mets);
    rxnDEgenes = [];
    for k=1:length(DEgUP)
            x = find(model.rxnGeneMat(:,find(strcmp(model.genes,DEgUP(k)))));
            rxnDEgenes = [rxnDEgenes;x];
    end
    rxnDEgenes = unique(rxnDEgenes);
    objUP = [];
    for j=1:length(metsUP)
        rxns   = find(model.S(b(j),:));
        indexes = intersect(rxns,rxnDEgenes);
        objUP = [objUP;indexes];
    end 
    objUP = unique(objUP);


    [a,b] = ismember(metsDN,model.mets);
    rxnDEgenes = [];
    for k=1:length(DEgDN)
            x = find(model.rxnGeneMat(:,find(strcmp(model.genes,DEgDN(k)))));
            rxnDEgenes = [rxnDEgenes;x];
    end
    rxnDEgenes = unique(rxnDEgenes);
    objDN = [];
    for j=1:length(metsDN)
        rxns   = find(model.S(b(j),:));
        indexes = intersect(rxns,rxnDEgenes);
        objDN = [objDN;indexes];
    end 
    objDN = unique(objDN);
    cVect = zeros(length(model.rxns),1);
    cUP = cVect;
    cUP(objUP) = 1;
    cDN = cVect;
    cDN(objDN) = -1;
    cVect = cUP+cDN;
    cIdxs = find(cVect);
    %load('../../models/candida_intermedia/cintGEM_curated.mat')

    temp = changeMedia_batch(model,'D-galactose exchange',1);
    temp = setParam(temp,'obj','r_4041',1);
    sol = solveLP(temp,1);
    temp = setParam(temp,'lb','r_4041',0.9*sol.x(find(temp.c)));
    temp = setParam(temp,'obj',1:length(model.rxns),cVect);
    sol = solveLP(temp,1);
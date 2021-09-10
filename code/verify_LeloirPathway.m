%explore lactose metbolism in C. intermedia model
yeastGEM = load('yeastGEM.mat');
yeastGEM = ravenCobraWrapper(yeastGEM.model);

load('Candida_intermedia.mat')
model = ravenCobraWrapper(reducedModel);
model.proteins = reducedModel.proteins;
%substitute gene IDs
model = substituteGENEids(model);
leloirGenes = {'YBR019C' 'YBR020W' 'YBR018C' 'YHL012W' 'YKL127W' 'YCL040W'};
leloirShort = {'GAL10' 'GAL1' 'GAL7' 'YHL012W' 'PGM1' 'GLK1'};

%Check the presence of leloir's genes in cintGEM, print rxn formulas
notPresent = [];
for i=1:length(leloirGenes)
    idx = find(strcmpi(model.genes,leloirGenes{i}));
    if ~isempty(idx)
        for j=1:length(idx)
            disp(leloirGenes{i})
            rxnIdxs = find(model.rxnGeneMat(:,idx(j)));
            printModel(model,rxnIdxs)
        end
    else
        disp(['Gene: ' leloirGenes{i} ' is not present in cIntGEM']) 
        notPresent = [notPresent;leloirGenes(i)];
    end
end
%It was found that reversibility for one of the reactions encoded by GAL7
%is wrong, let's fix it! (according to KEGG)
% YBR018C -> r_0459
idx = find(strcmpi(model.rxns,'r_0459'));
model.lb(idx) = -1000;
model.rev(idx) = 1;
%It was found that one of the reactions encoded by GAL10 is missing in Cint
%let's introduce it 
%obtain grRule for the GAL10 encoded rxns (already present)
grRule = model.grRules{find(strcmpi(model.rxns,'r_1070'))};
newRxns = {'D-galactose[c] => alpha-D-Galactose[c]' ...
          'lactose[e] => '};
rxnsToAdd.equations = newRxns; 
% Define reaction names
rxnsToAdd.rxns     = {'galMut' 'lac_ex'};
rxnsToAdd.rxnNames = {'galactose mutarotase' 'lactose exchange'};
% Define objective and bounds
rxnsToAdd.c  = [0 0];
rxnsToAdd.lb = [0 -1000];
rxnsToAdd.ub = [1000 0];
rxnsToAdd.grRules = {grRule ''};
model_leloir = addRxns(model,rxnsToAdd,3);
%It wa also found that conversion from lactose to D-galactose (r_5119) is defined in
%the reverse direction in this model
idx = find(strcmpi(model_leloir.rxns,'r_5119'));
printModel(model_leloir,idx)
model_leloir.rev(idx) = 1;
model_leloir.lb(idx) = -1000;
model_leloir.ub(idx) = 1000;
printModel(model_leloir,idx)

%now let's just see where the lac exchange is! Index 3988
lacEx = find(strcmpi(model_leloir.rxns,'lac_ex'));
printModel(model_leloir,lacEx)


model = changeMedia(model_leloir,1);


save('cint_leloir_AVR.mat','model_leloir')

 


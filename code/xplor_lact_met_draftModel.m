%explore lactose metbolism in C. intermedia model
load('../models/candida_intermedia/Candida_intermedia.mat')
model = ravenCobraWrapper(reducedModel);
model.proteins = reducedModel.proteins;
%Find lactose metabolites
lacIdxs = find(contains(model.metNames,'lactose'));
lacTable = table(model.metNames(lacIdxs),model.mets(lacIdxs), lacIdxs,...
                model.compNames(model.metComps(lacIdxs)),...
                'VariableNames',{'metNames' 'mets' 'idxs' 'compartments'});
%Find all related reactions for "lactose" entries
lacTable.rxnIdxs = cell(height(lacTable),1);
lacTable.rxnIds = cell(height(lacTable),1);
lacTable.rxnNames = cell(height(lacTable),1);
lacTable.associatedGenes = cell(height(lacTable),1);
for i=1:length(lacIdxs)
    rxnIdxs = find(model.S(lacIdxs(i),:));
    lacTable.rxnIdxs{i} = rxnIdxs;
    lacTable.rxnIds{i} = model.rxns(rxnIdxs);
    lacTable.rxnNames{i} = model.rxnNames(rxnIdxs);
    lacTable.associatedGenes{i} = model.grRules(rxnIdxs);
    %Translate all grRules here to cerevisiae IDs
end
%Exchange reactions table
[a,b] = getExchangeRxns(model);
exchTable = table(a,model.rxnNames(b),model.grRules(b),'VariableNames',{'rxns' 'rxnNames' 'grRules'});
exchTable.exch_mets = cell(height(exchTable),1);
exchTable.exch_metNames = cell(height(exchTable),1);
for i=1:length(b)
    metIdxs = find(model.S(:,b(i)));
    exchTable.exch_mets(i) = model.mets(metIdxs);
    exchTable.exch_metNames(i) = model.metNames(metIdxs);
end
x = find(contains(exchTable.exch_metNames,'lactose'));
%There's no exchange reaction for lactose in this model, as expected, let's
%introduce it!

% Define reactions equations
exchange_lac   = 'lactose[e] => ';
rxnsToAdd.equations = {exchange_lac}; 
% Define reaction names
rxnsToAdd.rxns     = {'lac_ex'};
rxnsToAdd.rxnNames = {'lactose exchange'};
% Define objective and bounds
rxnsToAdd.c  = [0];
rxnsToAdd.lb = [-1000];
rxnsToAdd.ub = [1000];

model_lac = addRxns(model,rxnsToAdd,3);
%Evaluate if rxn can carry flux
I = haveFlux(model_lac,1E-6,'lac_ex');
model = model_lac;
%the introduced reaction cannot carry flux, let's check the rest of the
%pathway
%The rxns of interest are reported in the lacTable,let's regenerate the
%table with the introduced reaction
lacIdxs = find(contains(model.metNames,'lactose'));
lacTable = table(model.metNames(lacIdxs),model.mets(lacIdxs), lacIdxs,...
                model.compNames(model.metComps(lacIdxs)),...
                'VariableNames',{'metNames' 'mets' 'idxs' 'compartments'});
%Find all related reactions for "lactose" entries
lacTable.rxnIdxs = cell(height(lacTable),1);
lacTable.rxnIds = cell(height(lacTable),1);
lacTable.rxnNames = cell(height(lacTable),1);
lacTable.associatedGenes = cell(height(lacTable),1);
lacTable.haveFlux = cell(height(lacTable),1);
for i=1:length(lacIdxs)
    rxnIdxs = find(model.S(lacIdxs(i),:));
    lacTable.rxnIdxs{i} = num2str(rxnIdxs);
    fluxVector = [];
    for j=1:length(rxnIdxs)
        ii = [num2str(haveFlux(model,1E-6,rxnIdxs(j))) ' / '];
        fluxVector = [fluxVector, ii];
    end
    lacTable.haveFlux{i} = fluxVector;
    lacTable.rxnIds{i} = strjoin(model.rxns(rxnIdxs),'/');
    lacTable.rxnNames{i} = strjoin(model.rxnNames(rxnIdxs),'/');
    lacTable.associatedGenes{i} = strjoin(model.grRules(rxnIdxs),'/');
end

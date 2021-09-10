%explore galactitol production in C. intermedia model
load('cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x)
%Model is already growing on lactose, growth is 0.20682
%We need to block the leloir genes first to see galactitol production!!


gtoh = find(contains(model.metNames,'galactitol'));
model.metNames(gtoh)
model.metComps(gtoh) %We already have gtoh(c) and gtoh(e)
sol.x(gtoh) %No flux in any of them
%Let's find the galactitol exchange reaction
gtoh_x = find(contains(model.rxnNames,'tol exchange'));
model.rxnNames(gtoh_x)
printModel(model,gtoh_x) %D-glucitol[e] =>  [0 1000]

lac_x = find(contains(model.rxnNames,'lactose exchange'));
model.rxnNames(lac_x) %lactose[e] <=>  [-1 1000]
printModel(model,lac_x)


%it looks like it doesn't exist in the model! Let's add it
newRxns = {'galactitol[c] <=> galactitol[e]'};
rxnsToAdd.equations = newRxns;
rxnsToAdd.rxnNames = {'galactitol exchange'};
rxnsToAdd.rxns     = {'gtoh_ex'};
%Define objective and bounds
rxnsToAdd.c  = 0;
rxnsToAdd.lb = 1000;
rxnsToAdd.ub = 1000;
% %genes to add
rxnsToAdd.grRules = {'xyl1_2'};
% Introduce changes to the model
model_oxido = addRxns(model_oxido,rxnsToAdd,3);



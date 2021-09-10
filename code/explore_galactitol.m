%explore galactitol production in C. intermedia model
load('cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x)
%Model is already growing on lactose, growth is 0.20682
%We need to block the leloir genes first to see galactitol production!!






%Let's find the galactitol exchange reaction
gtoh = find(contains(model.rxns,'ald_red_NADH'));
printModel(model,gtoh) %not exchange, but formation! we're getting close






gtoh = find(contains(model.rxnNames,'ol exchange'));
model.rxnNames(gtoh) 
%apparently it doesn't exist! S: Let's add it:
newRxns = {'D-galactose[c] + NADPH[c] => galactitol[c] + NADP(+)[c]'... 
           'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> L-sorbose[c] + NADP(+)[c]'};...
           %'D-Fructose[c] + ATP[c] <=> D-Fructose-6-Phosphate[c] + ADP[c]'};
rxnsToAdd.equations = newRxns;

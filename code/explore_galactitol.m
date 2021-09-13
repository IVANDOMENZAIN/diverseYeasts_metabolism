%explore galactitol production in C. intermedia model
load('cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x)
%Model is already growing on lactose, growth is 0.20682

%it looks like other yeasts (such as galactose-fermenting S.
%cerevisiae strains) display intracellular accumulation of galactitol (https://doi.org/10.1002/bit.21890), so

% Max growth on lactose = 27.0712
% Max growth on glucose = 27.0675 WHY. WHY. WHY. WHY.

lacEx = find(strcmpi(model.rxns,'lac_ex'));
printModel(model,lacEx)
model.lb(lacEx)   = 0;
model.ub(lacEx)   = 1000;

glu_ex = find(contains(model.rxnNames,'glucose exchange'));
printModel(model,glu_ex)
model.lb(glu_ex)   = -1000;
model.ub(glu_ex)   = 1000;

sol = solveLP(model,1);
printFluxes(model,sol.x)

%Let's see if galactitol is produced
gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
sol.x(gtoh_x)
%No galactitol produced at all

%Let's just see if we have flux along the redox pathway
rxns = {'ald_red_NADPH' 'r_4983' 'xyl_hex_red' 'r_5174' 'r_0323'};
fluxes=haveFlux(model,1E-6,rxns);
% They do

%let's see what happens when we block the leloir pathway:

index = find(contains(model.rxns,'r_0458')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 1000;
index = find(contains(model.rxns,'r_4222')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 1000;
sol = solveLP(model,1);
printFluxes(model,sol.x) %Growth: 0.19524

%Let's see if galactitol is produced now
gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
sol.x(gtoh_x)
%It was! Flux was 1

%what happens if we up the lactose consumption?
model.lb(3988)   = -1000;
model.ub(3988)   = 1000;
sol = solveLP(model,1);
printFluxes(model,sol.x)
sol.x(gtoh_x) %We have the same quantity of lactose being converted into galactitol, so maybe it's not being utilized further??
%The other reaction that utilizes galactitol is r_4983, let's see how it
%carries flux
gtoh_cons = find(contains(model.rxns,'r_4983'));
sol.x(gtoh_cons) %Same flux, the model is consuming all of the galactitol as soon as it's produced.
% We need to figure out if there is an exchange reaction that excretes the
% galactitol or if there's a bottleneck in r_4983 and it's being stuck inside of the cell

%Let's also run flux through the whole pathway
%The entire pathway just keeps moving forward!
OR_path = find(contains(model.rxnNames,'L-xylo-3-hexulose reductase'))
sol.x(OR_path)






































% gtoh = find(contains(model.metNames,'galactitol'));
% model.metNames(gtoh)
% model.metComps(gtoh) %We already have gtoh(c) and gtoh(e)
% sol.x(gtoh) %No flux in any of them
% %Let's find the galactitol exchange reaction
% gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
% model.rxnNames(gtoh_x)
% printModel(model,gtoh_x) %D-glucitol[e] =>  [0 1000]
% 
% lac_x = find(contains(model.rxnNames,'lactose exchange'));
% model.rxnNames(lac_x) %lactose[e] <=>  [-1 1000]
% printModel(model,lac_x)
% 
% % %it looks like it doesn't exist in the model! Let's add it
% % newRxns = {'galactitol[e] <=>'};
% % rxnsToAdd.equations = newRxns;
% % rxnsToAdd.rxnNames = {'galactitol exchange'};
% % rxnsToAdd.rxns     = {'gtoh_ex'};
% % %Define objective and bounds
% % rxnsToAdd.c  = 0;
% % rxnsToAdd.lb = 0;
% % rxnsToAdd.ub = 1000;
% % % %genes to add
% % rxnsToAdd.grRules = {'xyl1_2'};
% % % Introduce changes to the model
% % model_oxido = addRxns(model_oxido,rxnsToAdd,3);
% 
% %Before we start adding reactions, 

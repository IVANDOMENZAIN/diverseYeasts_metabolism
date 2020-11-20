%explore lactose metbolism in C. intermedia model
load('../models/candida_intermedia/cint_leloir.mat')

%Adding reactions from aspergillus niger oxidoreductive pathway. Reference:
%10.1016/j.fbr.2007.02.006
% Define reactions equations

%introducing the oxidoreductive pathway reactions from A.
%niger/A.nidulans/T.reesei. Reference:﻿﻿10.1074/jbc.M112.372755
newRxns = {'D-galactose[c] + NADPH[c] => galactitol[c] + NADP(+)[c]'... 
           'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> L-sorbose[c] + NADP(+)[c]'};...
           %'D-Fructose[c] + ATP[c] <=> D-Fructose-6-Phosphate[c] + ADP[c]'};
rxnsToAdd.equations = newRxns;

% Define reaction names
rxnsToAdd.rxnNames = {'aldose reductase' 'L-xylo-3-hexulose reductase'};% 'hexokinase'};
rxnsToAdd.rxns     = {'ald_red' 'xyl_hex_red'};% 'hxk'};
%Define objective and bounds
rxnsToAdd.c  = [0 0];
rxnsToAdd.lb = [0 -1000];
rxnsToAdd.ub = [1000 1000];
% %genes to add
genesToAdd.genes          = {'xyl1' 'G0RNA2'};
genesToAdd.geneShortNames = {'xyl1' 'lxr4'};  

rxnsToAdd.grRules         = {'xyl1' 'G0RNA2'};
%LEt's evaluate biomass production before integrating the pathway

model = changeMedia(model,1);
sol1 = solveLP(model,1);
printFluxes(model,sol1.x)
% Introduce changes to the model
model_oxido = addGenesRaven(model,genesToAdd);
model_oxido = addRxns(model_oxido,rxnsToAdd,3);
%Evaluate if rxn can carry flux
I = haveFlux(model_oxido,1E-6,'ald_red');
I2 = haveFlux(model_oxido,1E-6,'xyl_hex_red');
%LEt's evaluate biomass production
model = changeMedia(model_oxido,1);
sol2 = solveLP(model,1);
printFluxes(model,sol2.x)


%the introduced reaction cannot carry flux, let's check the rest of the
%pathway

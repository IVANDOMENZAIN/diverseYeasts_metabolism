%explore lactose metbolism in C. intermedia model
load('/Users/ramana/Documents/GitHub/diverseYeasts_metabolism/code/models/Candida_intermedia.mat')
model = ravenCobraWrapper(reducedModel);
model.proteins = reducedModel.proteins;
%Adding reactions from aspergillus niger oxidoreductive pathway. Reference:
%10.1016/j.fbr.2007.02.006
% Define reactions equations
%MCR_rxn        = 'malonyl-CoA[c] + 2 NADPH[c] => 3-hydroxypropionic acid[c] + 2 NADP(+)[c]';
%transport_3HP  = '3-hydroxypropionic acid[c] => 3-hydroxypropionic acid[e]';
%introducing the oxidoreductive pathway reactions from A.
%niger/A.nidulans/T.reesei. Reference:﻿﻿10.1074/jbc.M112.372755
newRxns = {'D-galactose[c] + NADPH[c] => galactitol[c] + NAD(+)[c]'... 
    'galactitol[c] + NAD(+)[c] <=> L-xylo-3-hexulose[c] + NADH[c]' ...
         'L-xylo-3-hexulose[c] + NADPH[c] <=> D-Sorbitol[c] + NADH[c]'...
         'D-Sorbitol[c] + NAD(+)[c] <=> D-Fructose[c] + NADH[c]'...
         'D-Fructose[c] + ATP[c] <=> D-Fructose-6-Phosphate[c] + ADP[c]';
rxnsToAdd.equations = newRxns;

% Define reaction names
rxnsToAdd.rxnNames     = {'aldose reductase' 'galactitol dehydrogenase' 'L-xylo-3-hexulose reductase' 'D-sorbitol dehydrogenase' 'hexokinase'};
rxnsToAdd.rxns         = {'ald_red' 'gal_deh' 'xyl_hex_red' 'sorb_deh' 'hxk'};
%Define objective and bounds
rxnsToAdd.c  = [0 0 0 0 0 0];
rxnsToAdd.lb = [-1000 -1000 -1000 -1000 -1000 -1000];
rxnsToAdd.ub = [1000 1000 1000 1000 1000 1000];

% % Metabolites to Add
metsToAdd.mets          = {'gal_tol' 'xylo_hex' 's_ose' 'S_tol''D-fructose'};
%randomly named the metabolites
metsToAdd.metNames      = {'galctitol' 'D-Sorbose' 'D-Sorbitol' 'D-Fructose-6-Phosphate' 'L-xylo-6-hexulose'};
% metsToAdd.compartments  = {'c' 'e'};
% %genes to add
genesToAdd.genes          = {'xyl1' 'ladA' 'xhrA' 'sdhA' 'hxk'};
% genesToAdd.geneShortNames = {'MCR'};
% not really sure about what I had to add here.  
% rxnsToAdd.grRules         = {'MCR' '' ''};
% Introduce changes to the model
model_oxido = addGenesRaven(model,genesToAdd);
model_oxido = addMets(model_oxido,metsToAdd);
model_oxido = addRxns(model,rxnsToAdd,3);
%Evaluate if rxn can carry flux
I = haveFlux(model_oxido,1E-6,'');
model = model_oxido;
%the introduced reaction cannot carry flux, let's check the rest of the
%pathway

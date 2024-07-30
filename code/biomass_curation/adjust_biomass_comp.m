% Kamesh Peri.      Last update: 2024-03-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../../models/candida_intermedia/cintGEM_oxido.mat')
%%correct reaction
% * L-xylo-3-hexulose reductase: 'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> L-sorbose[c] + NADP(+)[c]'};... G0RNA2 (lxr4)
%by this:  * L-xylo-3-hexulose reductase: 'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> D-glucitol[c] + NADP(+)[c]'};... G0RNA2 (lxr4)
rxnName = 'xyl_hex_red';
x = find(strcmp(model.rxns,'xyl_hex_red'));
products = find(model.S(:,x)>0);
disp(model.metNames(products))
products = products(1);
disp(model.metNames(products))
model.S(products,x) = 0;
y = find(strcmpi(model.metNames,'D-glucitol')); 
model.metComps(y)
y = y(1);
model.S(y,x) = 1;
constructEquations(model,x,true)


rxnsToAdd = struct();
newRxns = {'galactitol[c] + NAD[c] => D-tagatose[c] + NADH[c] + H+[c]'};...
           %'D-Fructose[c] + ATP[c] <=> D-Fructose-6-Phosphate[c] + ADP[c]'};
rxnsToAdd.rxns = {'galactitol_dhd'};
rxnsToAdd.equations = newRxns;

% Define reaction names
rxnsToAdd.rxnNames = {'galactitol dehydrogenase'};% 'hxk'};
%Define objective and bounds
rxnsToAdd.c  = [0];
rxnsToAdd.lb = [0];
rxnsToAdd.ub = [1000];
% %genes to add
genesToAdd.genes          = {'galactitol_dh'};
genesToAdd.geneShortNames = {'galactitol_dh'};  
rxnsToAdd.grRules         = {'galactitol_dh'};
% Introduce changes to the model
model_oxido = addGenesRaven(model,genesToAdd);
model_oxido.proteins = [model_oxido.proteins; genesToAdd.genes'];
model_oxido = addRxns(model_oxido,rxnsToAdd,3);
%the model does not have any reaction for secreting galactitol, introduce
%it
newRxns = {'galactitol[e] <=> '};
rxnsToAdd.equations = newRxns; 
% Define reaction names
rxnsToAdd.rxns     = {'galactitol exchange'};
rxnsToAdd.rxnNames = {'galactitol exchange'};
% Define objective and bounds
rxnsToAdd.c  = [0];
rxnsToAdd.lb = [0];
rxnsToAdd.ub = [1000];
rxnsToAdd.grRules = {''};
model_oxido = addRxns(model_oxido,rxnsToAdd,3);
model = model_oxido;

%correct reveersibilities of lactose metabolism, start with leloir
%[~,b] = ismember(galGenes,model.genes);
model = setParam(model,'lb','r_0459',-1000); %GAL7
model = setParam(model,'lb','galMut',-1000); %GAL10 (mutarotase part)
model.rev(find(strcmp(model.rxns,'galMut'))) = 1;
%now for Ox-red
x = find(strcmp(model.rxns,'ald_red_NADH'));
model.rev(x) =0;
model.lb(x) = 0;
x = find(strcmp(model.rxns,'ald_red_NADPH'));
model.rev(x) =0;
model.lb(x) = 0;

%verify growth on lactose
disp('Growth on glucose')
model = changeMedia_batch(model,'D-glucose exchange',1);
sol = solveLP(model,1);
printFluxes(model,sol.x,true)
disp(' ')
disp('Growth on lactose')
model = changeMedia_batch(model,'lactose exchange',1);
sol = solveLP(model,1);
printFluxes(model,sol.x,true)
disp(' ')
disp('Growth on D-galactose')
model = changeMedia_batch(model,'D-galactose exchange',1);
sol = solveLP(model,1);
printFluxes(model,sol.x,true)
disp(' ')
disp('Growth on lactose')
model = changeMedia_batch(model,'lactose exchange',1);
model = setParam(model,'obj',3736,1); %GAL7
sol = solveLP(model,1);
printFluxes(model,sol.x,true)
disp(' ')
pause
%from chemostat data it was found that the GUR at 0 dilution rate must
%correspond to 0.03 mmol/gDw h, fix this GUR and max. NGAM to obtain its LB
x = find(strcmpi(model.rxnNames,'non-growth associated maintenance reaction'));
model = changeMedia_batch(model,'D-glucose exchange',0.03);
model.lb(x) = 0;
model.ub(x) = 1000;
temp = setParam(model,'obj',x,1);
sol = solveLP(temp);
model.lb(x) = sol.x(x);
%correct stoichiometry in complex I, lets start with the theoretical valuer of
% 2 (as a basis coeff. for proton translocation)
%introducing the oxidoreductive pathway reactions from A.
%niger/A.nidulans/T.reesei. Reference:﻿﻿10.1074/jbc.M112.372755
model = changePOratio(model,2.2);

Ptot = 0.438; %average across chemostats in g protein / gCDW

posBiomass = find(contains(model.rxnNames,'biomass'));
constructEquations(model,posBiomass)
%the biomass rxn is the modular one, print the pseudoreaction for each of
%the different modular components (carbs, lipids, prots, etc.
% posProt = find(contains(model.rxnNames,'rotein pseudoreaction'));
% constructEquations(model,posProt)
% posLip = find(contains(model.rxnNames,'ipid backbone pseudoreaction'));
% constructEquations(model,posLip)
% posCarb = find(contains(model.rxnNames,'arbohydrate pseudoreaction'));
% constructEquations(model,posCarb)
% posIon = find(contains(model.rxnNames,'ion pseudoreaction'));
% constructEquations(model,posIon)
% posRNA = find(contains(model.rxnNames,'RNA pseudoreaction'));
% constructEquations(model,posRNA)
% posDNA = find(contains(model.rxnNames,'DNA pseudoreaction'));
% constructEquations(model,posDNA)

%Components of biomass:  (from yeastGEM)
%        id         MW [g/mol]  class     name
comps = {'s_0404'	89.09       'P'     % A     Alanine         ala
         's_0542'	121.16      'P'     % C     Cysteine        cys
         's_0432'   133.11      'P'     % D     Aspartic acid   asp
         's_0748'   147.13      'P'     % E     Glutamic acid   glu
         's_1314'   165.19      'P'     % F     Phenylalanine   phe
         's_0757'   75.07       'P'     % G     Glycine         gly
         's_0832'   155.15      'P'     % H     Histidine       his
         's_0847'   131.17      'P'     % I     Isoleucine      ile
         's_1099'   146.19      'P'     % K     Lysine          lys
         's_1077'   131.17      'P'     % L     Leucine         leu
         's_1148'   149.21      'P'     % M     Methionine      met
         's_0430'   132.12      'P'     % N     Asparagine      asn
         's_1379'   115.13      'P'     % P     Proline         pro
         's_0747'   146.14      'P'     % Q     Glutamine       gln
         's_0428'   174.2       'P'     % R     Arginine        arg
         's_1428'   105.09      'P'     % S     Serine          ser
         's_1491'   119.12      'P'     % T     Threonine       thr
         's_1561'   117.15      'P'     % V     Valine          val
         's_1527'   204.23      'P'     % W     Tryptophan      trp
         's_1533'   181.19      'P'     % Y     Tyrosine        tyr
         's_0001'	180.16      'C'     % (1->3)-beta-D-glucan
         's_0004'	180.16      'C'     % (1->6)-beta-D-glucan
         's_0773'   180.16      'C'     % glycogen
         's_1107'   180.16      'C'     % mannan
         's_1520'   342.296 	'C'     % trehalose
         's_0423'   347.22      'R'     % AMP
         's_0526'   323.2       'R'     % CMP
         's_0782'   363.22      'R'     % GMP
         's_1545'   324.18      'R'     % UMP
         's_0584'   331.22      'D'     % dAMP
         's_0589'   307.2       'D'     % dCMP
         's_0615'   345.21      'D'     % dGMP
         's_0649'   322.21      'D'     % dTMP
         's_3714'   852.83      'N'     % heme a
         's_1405'   376.36      'N'     % riboflavin
         's_1467'   96.06       'N'};   % sulphate

%WITH THIS we have checked that all the metabolites composing the biomass
%reaction are all consistent with those in yeastGEM

%Get main fractions:
[Pbase,X] = getFraction(model,comps,'P',0);
[Cbase,X] = getFraction(model,comps,'C',X);
[Rbase,X] = getFraction(model,comps,'R',X);
[Dbase,X] = getFraction(model,comps,'D',X);
[Lbase,X] = getFraction(model,comps,'L',X);
%The sum Protein + carbohydrates + lipid backbones + RNA + DNA 
%accounts for 99.75% of the biomass in the model

%Now rescale biomass according to the provided experimental value of total
%protein content. This will rescale carbs and lipids using the same
%proportions as in the base biomass equation
Ctot = Cbase + (Pbase - Ptot)*Cbase/(Lbase+Cbase);
Ltot = Lbase + (Pbase - Ptot)*Lbase/(Lbase+Cbase);
%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fL = Ltot/Lbase;
%Change compositions:
modelMod = rescalePseudoReaction(model,'protein',fP);
modelMod = rescalePseudoReaction(modelMod,'carbohydrate',fC);
%If model contain SLIMER reactions (separate pseudoreactions for
%lipid chains and backbones
modelMod = rescalePseudoReaction(modelMod,'lipid backbone',fL);
modelMod = rescalePseudoReaction(modelMod,'lipid chain',fL);
%Check how stoichiometries have changed for each of the biomass components
% constructEquations(modelMod,posBiomass)
% constructEquations(model,posProt)
% constructEquations(modelMod,posProt)
% constructEquations(model,posCarb)
% constructEquations(modelMod,posCarb)
% constructEquations(model,posLip)
% constructEquations(modelMod,posLip)
%recompute the sum of mass fractions (Protein + carbohydrates + lipid backbones + RNA + DNA) 
[~,X] = getFraction(modelMod,comps,'P',0);
[~,X] = getFraction(modelMod,comps,'C',X);
[~,X] = getFraction(modelMod,comps,'R',X);
[~,X] = getFraction(modelMod,comps,'D',X);
[~,X] = getFraction(modelMod,comps,'L',X);
clc

for j=1:1
    GAM = fitGAM(modelMod);
    modelMod =changeGAM(modelMod,GAM);
    POratio  = fitPOratio(modelMod);
    modelMod = changePOratio(modelMod,POratio);
end

%save curated model
model = modelMod;
model = changeMedia_batch(model,'lactose exchange',1);
model = setParam(model,'obj',3736,1);
%block ATP:D-tagatose 6-phosphotransferase
model = setParam(model,'ub','r_4393',0);
model = setParam(model,'lb','r_4393',0);


sol = solveLP(model,1);
printFluxes(model,sol.x,true)
save('../../models/candida_intermedia/cintGEM_curated.mat','model')


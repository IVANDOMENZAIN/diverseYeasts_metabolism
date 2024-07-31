function output = adjust_biomass_comp
clear
clc
load('../../models/candida_intermedia/cintGEM_gene_curated.mat')
%%correct reaction
% * L-xylo-3-hexulose reductase: 'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> L-sorbose[c] + NADP(+)[c]'};... G0RNA2 (lxr4)
%by this:  * L-xylo-3-hexulose reductase: 'L-xylo-3-hexulose[c] + NADPH[c] + H+[c] <=> D-glucitol[c] + NADP(+)[c]'};... G0RNA2 (lxr4)
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

%add galactitol dehydrogenase (based on evidence on wild-type growth on
%galactose
rxnsToAdd = struct();
newRxns = {'galactitol[c] + NAD[c] => D-tagatose[c] + NADH[c] + H+[c]'};...
rxnsToAdd.rxns = {'galactitol_dhd'};
rxnsToAdd.equations = newRxns;
% Define reaction names
rxnsToAdd.rxnNames = {'galactitol dehydrogenase'};
%Define objective and bounds
rxnsToAdd.c  = 0;
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;
% %genes to add
rxnsToAdd.grRules = {'Seq_2272'};
% Introduce changes to the model
model_oxido = model;%addGenesRaven(model,genesToAdd);
model_oxido = addRxns(model_oxido,rxnsToAdd,3);
%the model does not have any reaction for secreting galactitol, introduce
%it (based on secretion/uptake experimental phenotype
newRxns = {'galactitol[e] <=> '};
rxnsToAdd.equations = newRxns; 
% Define reaction names
rxnsToAdd.rxns     = {'galactitol exchange'};
rxnsToAdd.rxnNames = {'galactitol exchange'};
% Define objective and bounds
rxnsToAdd.c  = 0;
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;
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

%from chemostat data it was found that the GUR at 0 dilution rate must
%correspond to 0.03 mmol/gDw h, fix this GUR and max. NGAM to obtain its LB
x = find(strcmpi(model.rxnNames,'non-growth associated maintenance reaction'));
model = changeMedia_batch(model,'D-glucose exchange',0.03);
model.lb(x) = 0;
model.ub(x) = 1000;
temp = setParam(model,'obj',x,1);
sol = solveLP(temp);
model.lb(x) = sol.x(x);

%calibrate biomass comp
modelMod = calibrate_biomass_pseudoreaction(model);

%correct stoichiometry in complex I, lets start with the theoretical valuer of
% 2.5 (as a basis coeff. for proton translocation, and in agreement with
% theory for organisms with complex I in ETC)
Theor_PO = 2.5;
modelMod = changePOratio(modelMod,Theor_PO);
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
%model = setParam(model,'ub','r_4393',0);
%model = setParam(model,'lb','r_4393',0);
%
sol = solveLP(model,1);
printFluxes(model,sol.x,true)
save('../../models/candida_intermedia/cintGEM_curated.mat','model')
output = model;
%generate version-controllable files
formulas = constructEquations(model);
rxns = model.rxns;
rxnNames = model.rxnNames;
grRules = model.grRules;
modelTable = table(rxns,rxnNames,formulas, grRules);
writetable(modelTable,'../../models/candida_intermedia/cintGEM_curated.txt','WriteVariableNames',true,'Delimiter','\t','QuoteStrings',false);

%add version control
genes = model.genes;
shortnames = model.geneShortNames;
orthologues = model.orthologues;
proteins = model.proteins;
gene_table = table(genes,shortnames,orthologues,proteins);
writetable(gene_table,'../../models/candida_intermedia/gene_table_curated.txt','Delimiter','\t','QuoteStrings',false);
end

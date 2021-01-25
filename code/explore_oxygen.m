%Explore Oxygen usage by model
%C. intermedia is a typical aerobic
%(https://www.lgcstandards-atcc.org/products/all/14439.aspx?geo_country=se#culturemethod)
%We load the OR model
load('../models/candida_intermedia/cintGEM_oxido.mat');
%Find O2 exchange
O2_ex = find(strcmpi(model.rxnNames,'oxygen exchange'))
O2_ex = model.rxns(O2_ex);
%Block oxygen exchange
model = setParam(model,'lb',O2_ex,0);
model = setParam(model,'ub',O2_ex,0);
%Run FBA and print fluxes
sol = solveLP(model,1)
printFluxes(model,sol.x)
%No growth! 

%Trying with Leloir model
load('../models/candida_intermedia/cint_leloir.mat');
%Find O2 exchange
O2_ex = find(strcmpi(model.rxnNames,'oxygen exchange'))
O2_ex = model.rxns(O2_ex);
%Block oxygen exchange
model = setParam(model,'lb',O2_ex,0);
model = setParam(model,'ub',O2_ex,0);
%Run FBA and print fluxes
sol2 = solveLP(model,1)
printFluxes(model,sol2.x)
%Still no growth!

%Checking model for anaerobic growth on YeastGEM (model =
%anaerobicModel(model))
%First checking if we have the necessary metabolites and exch reactions
find(strcmp(model.mets,'s_3714[c]')) 
find(strcmp(model.mets,'s_1198[c]'))
find(strcmp(model.mets,'s_1203[c]'))
find(strcmp(model.mets,'s_1207[c]'))
find(strcmp(model.mets,'s_1212[c]'))
find(strcmp(model.mets,'s_0529[c]'))
%Didn't find anything, will look without the compartment suffix
find(strcmp(model.mets,'s_3714')) 
find(strcmp(model.mets,'s_1198'))
find(strcmp(model.mets,'s_1203'))
find(strcmp(model.mets,'s_1207'))
find(strcmp(model.mets,'s_1212'))
find(strcmp(model.mets,'s_0529'))
%Got them!
find(strcmp(model.rxns,'r_1994')) %palmitoleate
find(strcmp(model.rxns,'r_2106')) %zymosterol
find(strcmp(model.rxns,'r_2134')) %14-demethyllanosterol
find(strcmp(model.rxns,'r_2137')) %ergosta-5,7,22,24(28)-tetraen-3beta-ol
find(strcmp(model.rxns,'r_2189')) %oleate
find(strcmp(model.rxns,'r_0713'))
find(strcmp(model.rxns,'r_0714')) %<--- couldn't find this one!
find(strcmp(model.rxns,'r_0487'))
find(strcmp(model.rxns,'r_0472'))
find(strcmp(model.rxnNames,'biomass pseudoreaction'))
%We have almost all of them!
find(strcmp(model.metNames,'ATP [cytoplasm]'))


model = anaerobicModel(model)

%No growth observed, running FSEOF if growth relies on oxygen or if we need
%to add something to the media hehe



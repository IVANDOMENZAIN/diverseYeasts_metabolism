%Explore Oxygen usage by model
%C. intermedia is a typical aerobic
%(https://www.lgcstandards-atcc.org/products/all/14439.aspx?geo_country=se#culturemethod)
%We load the OR model
load('../models/candida_intermedia/cintGEM_oxido.mat');
%Find O2 exchange
O2_ex = find(strcmpi(model.rxnNames,'oxygen exchange'))
O2_ex = model.rxns(O2_ex);
%Let's set minimal media conditions
model = changeMedia(model,'r_1714',1)
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
toRemove = {'s_3714' 's_1198' 's_1203' 's_1207' 's_1212' 's_0529'};
tempModel = model;
%Get index for biomass-related reactions
bioRxn = find(strcmpi(model.rxnNames,'biomass pseudoreaction'));
cofRxn = find(strcmpi(model.rxnNames,'cofactor pseudoreaction'));
proRxn = find(strcmpi(model.rxnNames,'protein pseudoreaction'));
ionRxn = find(strcmpi(model.rxnNames,'ion pseudoreaction'));
lipRxn = find(contains(model.rxnNames,'lipid chain'));
printModel(tempModel,bioRxn)
printModel(tempModel,cofRxn)
printModel(tempModel,proRxn)
printModel(tempModel,ionRxn)
printModel(tempModel,lipRxn)
%model.metNames(find(strcmp(model.mets,'s_3714[c]'))) 

printModel(model,find(strcmp(model.rxns,'r_1994'))) %palmitoleate
printModel(model,find(strcmp(model.rxns,'r_2106'))) %zymosterol
printModel(model,find(strcmp(model.rxns,'r_2134'))) %14-demethyllanosterol
printModel(model,find(strcmp(model.rxns,'r_2137'))) %ergosta-5,7,22,24(28)-tetraen-3beta-ol
printModel(model,find(strcmp(model.rxns,'r_2189'))) %oleate
printModel(model,find(strcmp(model.rxns,'r_0713')))
%printModel(find(strcmp(model.rxns,'r_0714'))) %<--- couldn't find this one!
printModel(model,find(strcmp(model.rxns,'r_0487')))
printModel(model,find(strcmp(model.rxns,'r_0472')))

toAdd = {'r_1994' 'r_2106' 'r_2134' 'r_2137' 'r_2189'};

for rxn=toAdd
    tempModel = setParam(tempModel,'lb',rxn,-1000);
    printModel(tempModel,find(strcmp(tempModel.rxns,rxn))) %palmitoleate
end
sol = solveLP(tempModel)
printFluxes(tempModel,sol.x)

toAdd = {'r_0713' 'r_0487' 'r_0472'};

for rxn=toAdd
    tempModel = setParam(tempModel,'lb',rxn,0);
    tempModel = setParam(tempModel,'ub',rxn,0);
    printModel(tempModel,find(strcmp(tempModel.rxns,rxn))) %palmitoleate
end
sol = solveLP(tempModel)
printFluxes(tempModel,sol.x)

%It seems that the Cint model cannot grow without oxygen (using anaerobic
%media composition for S. cerevisiae), so let's open all the possible and
%see if this is still the case
[a,b]= getExchangeRxns(tempModel);
exchNames = tempModel.rxnNames(b);
tempModel.lb(b) = -1000;
%close oxygen
%Block oxygen exchange
tempModel = setParam(tempModel,'lb',O2_ex,0);
tempModel = setParam(tempModel,'ub',O2_ex,0);

sol = solveLP(tempModel)
printFluxes(tempModel,sol.x)
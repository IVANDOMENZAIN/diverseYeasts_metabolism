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
%Set minimal media with glucose as carbon source
model = changeMedia(model,'r_1714',1);
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
exchRxns = tempModel.rxns(b);
tempModel.lb(b) = -1000;
%close oxygen
%Block oxygen exchange
tempModel = setParam(tempModel,'lb',O2_ex,0);
tempModel = setParam(tempModel,'ub',O2_ex,0);

sol = solveLP(tempModel)
printFluxes(tempModel,sol.x)
%There seems to be anaerobic growth if we open all of the exchange reactions
%Now, we will close each reaction separately to see which affect growth rate.

%toBlock = [];
%for i = 1:length(sol.x)
 %   if sol.x(i) < 0
  %      toBlock(i) = sol.x(i);
   % else
        
   % end
%end
%toBlock
for i = 1:length(exchRxns)
    %if sol.x(i)<0
     %   exchRxns(i)
    tempModel = setParam(tempModel,'lb',exchRxns(i),0);
    tempModel = setParam(tempModel,'ub',exchRxns(i),1000);
    tempModel = setParam(tempModel,'lb',O2_ex,0);
    tempModel = setParam(tempModel,'ub',O2_ex,0);
    sol = solveLP(tempModel);
    if abs(sol.f)<=1E-6
        disp('Essential component')
        fprintf("Growth without rxn = %g\n",sol.f);
        disp(exchNames(i))
    end
    tempModel = setParam(tempModel,'lb',exchRxns(i),-1000);
    tempModel = setParam(tempModel,'ub',exchRxns(i),1000);
  %  end
end

%We can try choosing 2 values at random
r = rand(length(exchRxns))

msize = numel(exchRxns);
exchRxns(randperm(msize, 2))
i = 0;
while i<=1000 
    tempModel = setParam(tempModel,'lb',exchRxns,-1000);
    x = exchRxns(randperm(msize, 2));
    tempModel = setParam(tempModel,'lb',x(1),0);
     tempModel = setParam(tempModel,'lb',x(2),0);
    tempModel = setParam(tempModel,'ub',x(2),1000);
    tempModel = setParam(tempModel,'ub',x(1),1000);
    tempModel = setParam(tempModel,'lb',O2_ex,0);
    tempModel = setParam(tempModel,'ub',O2_ex,0);
    sol = solveLP(tempModel);
    if abs(sol.f)<=1E-6
        disp('Essential component')
        fprintf("Growth without rxn = %g\n",sol.f);
        index = find(strcmp(exchRxns,x(1)));
        disp(exchNames(index))
        index = find(strcmp(exchRxns,x(2)));
        disp(exchNames(index))
    end
    %print(exchRxns(i), sol.f)
    tempModel = setParam(tempModel,'lb',x(1),-1000);
    tempModel = setParam(tempModel,'lb',x(2),-1000);
    i =i+1;
end
%Hard to get any conclusions

%Let's block component by component (in a cumulative way) and save those
%exchanges that when blocked tje growth rate is reduced
tempModel = setParam(tempModel,'lb',exchRxns,-1000);
tempModel = setParam(tempModel,'ub',exchRxns,1000);
tempModel = setParam(tempModel,'lb',O2_ex,0);
tempModel = setParam(tempModel,'ub',O2_ex,0);
model2    = tempModel;
sol = solveLP(tempModel);
growth = abs(sol.f);
essentialExchs = [];
allGrowthRates = [];
essential = {'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020'; ...  % potassium exchange
                    'r_4062'; ...
                    'r_4064'};  
[~,iB] = ismember(essential,exchRxns);
%exchRxns(iB) = {};
for i=1:length(exchRxns)
    if ~ismember(exchRxns{i},essential)
    tempModel = setParam(tempModel,'lb',exchRxns(i),0);
    sol = solveLP(tempModel);
    if abs(sol.f)<growth
        fprintf("Growth without rxn = %g\n",sol.f);
        disp(exchNames(i))
    end
    end
        growth = abs(sol.f);
    if ~isempty(growth)
        allGrowthRates = [allGrowthRates; growth];
    else
        allGrowthRates = [allGrowthRates; NaN];
    end
end
t = table(exchRxns,exchNames,allGrowthRates);
writetable(t,'../results/growthRates_exchRxns.txt','delimiter','\t')

%exploreCintPotential
load('../models/candida_intermedia/cintGEM_oxido.mat')
carbon_sources = {'D-glucose exchange' 'lactose exchange'};
growth_idx = find(strcmpi(model.rxnNames,'growth'));
%Test Model 
sol = solveLP(model,1);
printFluxes(model,sol.x)
gRate_ec = sol.x(growth_idx);
%Identify all metabolites that can be produced by the model
[a,b] = getExchangeRxns(model);
formulas = constructEquations(model,b);
production = table(a,model.rxnNames(b),b,formulas,'VariableNames',{'rxns' 'rxnNames' 'index' 'formulas'});
%different levels of biomass formation (relative to max growth)
bio_levels = [0 0.25 0.50 0.75 0.99];
%Iterate through each carbon source
for cSource = carbon_sources
    disp(cSource)
    %create a table of results including all information about exchange
    %reactions
    production = table(a,model.rxnNames(b),b,formulas,'VariableNames',{'rxns' 'rxnNames' 'index' 'formulas'});
    %Set minimal mineral media with a unit uptake rate of the desired
    %carbon source
    cS_id = model.rxns(find(strcmpi(model.rxnNames,cSource)));
    temp  = minimal_Y6(model,cS_id);
    %Get maximum growth rate subject to a unit carbon source uptake rate
    sol       = solveLP(temp,1);
    maxGrowth = sol.x(growth_idx);
    %for each suboptimal biomass level explore production potential for all
    %exchangeable metabolites
    for i=1:length(bio_levels)
        %set a lower bound of biomass formation subject to predefined
        %suboptimal level
        temp1 = setParam(temp,'lb',growth_idx,bio_levels(i)*maxGrowth);
        vector =zeros(height(production),1);
        %Iterate through all of the possible exchangeable metabolites
        for j=1:length(b)
            %Change objective function to the corresponding exchange
            %reaction
        	temp1 = setParam(temp1,'obj',b(j),1);
            %Unconstrain exchange rate
            temp1 = setParam(temp1,'ub',b(j),1000);
            %Solve!
            sol   = solveLP(temp1);
            vector(j) = sol.x(b(j));
        end
        %append vector of results to production table
        eval(['production.growth' num2str(i) '=vector;'])
    end
    %for each carbon source write results in .txt file
    if strcmpi(cSource,'D-glucose exchange')
        filename = '../results/Cint_glc_potential.txt';
    else
        filename = '../results/Cint_lac_potential.txt';
    end
    writetable(production,filename,'delimiter','\t','QuoteStrings',false)
end

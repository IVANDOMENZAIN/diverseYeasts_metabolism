function mod_data = simulateChemostat(model,exp_data,GAM)
%Modify GAM withouth changing the protein content:
parameters.exch_names{1} = 'growth';
parameters.exch_names{2} = 'D-glucose exchange';
parameters.exch_names{3} = 'oxygen exchange';
parameters.exch_names{4} = 'carbon dioxide exchange';
%Relevant positions:
exch_names  = parameters.exch_names;
pos(1)      = find(strcmp(model.rxnNames,exch_names{1}));
pos(2)      = find(strcmp(model.rxnNames,exch_names{2}));
pos(3)      = find(strcmp(model.rxnNames,exch_names{3}));
pos(4)      = find(strcmp(model.rxnNames,exch_names{4}));

%Simulate chemostats:
mod_data = zeros(size(exp_data));
for i = 1:length(exp_data(:,1))
    %Fix biomass
    exp_data(i,1)
    model = setParam(model,'lb',model.rxns(pos(1)),0.99*exp_data(i,1));
    %set an arbitrarily high glucose uptake rate
    model = setParam(model,'lb',model.rxns(pos(2)),-10);
    %provide oxygen
    model = setParam(model,'lb',model.rxns(pos(3)),-1000);
    %block CO2 uptake and allow its production
    model = setParam(model,'ub',model.rxns(pos(4)),1000);
    model = setParam(model,'lb',model.rxns(pos(4)),0);
    %minimize glucose
    model = setParam(model,'obj',model.rxns(pos(2)),1);
    sol   = solveLP(model,1);
    printFluxes(model,sol.x,true)
    pause
    %Store relevant variables:
    mod_data(i,:) = sol.x(pos)';
end
end
current = pwd;
%load models
model = importModel('../models/blastobotrys_mokoenaii.xml');
load('../models/ecBmo/ecModel_batch.mat')
%Test ecModel 
sol = solveLP(ecModel_batch,1);
growthIdx = find(ecModel_batch.c);
gRate_ec = sol.x(growthIdx);
%Rescale protein pool upper bound
newBound = (0.5/gRate_ec)*ecModel_batch.ub(end);
ecModel_batch.ub(end)  = newBound;
obj_ec = find(ecModel_batch.c);
sol = solveLP(ecModel_batch,1);
printFluxes(ecModel_batch, sol.x, true, 10^-6);
%Retrieve all exchange reactions
model = minimal_Y6(model,'r_1714');
solution = solveLP(model);
obj = (find(model.c));
GURidx = (find(strcmpi(model.rxns,'r_1714')));
gRate = solution.x(find(model.c));
ecModel_batch = setParam(ecModel_batch,'ub','r_1714_REV',1);
ecGURidx = (find(strcmpi(ecModel_batch.rxns,'r_1714_REV')));

sol = solveLP(ecModel_batch,1);
gRate_ec = sol.x(growthIdx);
[modelRxns,modelIdxs] = getExchangeRxns(model);
%Initialize variables for results
excRxnIds    = [];
modelRates   = [];
ecModelRates = [];
model_subRates   = [];
ecModel_subRates = [];
model_maxRates   = [];
ecModel_maxRates = [];
for i=1:length(modelRxns)
    ecM_idx = find(strcmpi(ecModel_batch.rxns,modelRxns{i}));
    if ~isempty(ecM_idx)
        tempModel   = setParam(model,'obj',modelIdxs(i),1);
        %Optimize for mximum production levels
        tempModel   = setParam(tempModel,'ub',modelIdxs(i),1000);
        tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
        solution    = solveLP(tempModel);
        solution_ec = solveLP(tempecModel);
        if ~isempty(solution.x) & ~isempty(solution_ec.x)
            %Optimize for growth and production
            tempModel   = setParam(tempModel,'lb',obj,0.5*gRate);
            tempecModel = setParam(tempecModel,'lb',obj_ec,0.5*gRate_ec);
            solution_sub    = solveLP(tempModel);
            solution_sub_ec = solveLP(tempecModel);
            if ~isempty(solution_sub.x) & ~isempty(solution_sub_ec.x)
                %Release GUR and optimize
                tempModel   = setParam(tempModel,'lb',obj,0.1);
                tempecModel = setParam(tempecModel,'lb',obj_ec,0.1);
                tempModel   = setParam(tempModel,'lb','r_1714',-1000);
                tempecModel = setParam(tempecModel,'ub','r_1714_REV',1000);
                solution_max   = solveLP(tempModel,1);
                solution_max_ec = solveLP(tempecModel,1);
                if ~isempty(solution_max.x) & ~isempty(solution_max_ec.x)
                    disp(['#' num2str(i) ' ' model.rxnNames{modelIdxs(i)}])
                    yield    = solution_max.x(modelIdxs(i))/abs(solution_max.x(GURidx));
                    yield_EC = solution_max_ec.x(ecM_idx)/abs(solution_max_ec.x(ecGURidx));
                    excRxnIds    = [excRxnIds; model.rxnNames(modelIdxs(i))];
                    modelRates   = [modelRates;solution.x(modelIdxs(i))];
                    ecModelRates = [ecModelRates;solution_ec.x(ecM_idx)];
                    model_subRates   = [model_subRates;solution_sub.x(modelIdxs(i))];
                    ecModel_subRates = [ecModel_subRates;solution_sub_ec.x(ecM_idx)];
                    model_maxRates   = [model_maxRates;yield];
                    ecModel_maxRates = [ecModel_maxRates;yield_EC];
                end
            end
        end
    end
end
results = table(excRxnIds,modelRates,ecModelRates,model_subRates,ecModel_subRates,model_maxRates,ecModel_maxRates);     
    

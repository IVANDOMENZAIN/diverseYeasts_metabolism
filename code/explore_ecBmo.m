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
excRxnIds = [];
rates_01  = [];
rates_25  = [];
rates_50  = [];
rates_75  = [];
rates_99  = [];

for i=1:length(modelRxns)
    ecM_idx = find(strcmpi(ecModel_batch.rxns,modelRxns{i}));
    if ~isempty(ecM_idx)
        tempModel   = setParam(model,'obj',modelIdxs(i),1);
        %Optimize for maximum production levels
        tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
        tempecModel = setParam(tempecModel,'lb',obj_ec,0.01*gRate_ec);
        solution_01 = solveLP(tempecModel);
        if ~isempty(solution_01.x)
            %Optimize for growth and production
            tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
            tempecModel = setParam(tempecModel,'lb',obj_ec,0.25*gRate_ec);
            solution_25 = solveLP(tempecModel);
            if ~isempty(solution_25.x)
                tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
                tempecModel = setParam(tempecModel,'lb',obj_ec,0.5*gRate_ec);
                solution_50 = solveLP(tempecModel);
                if ~isempty(solution_50.x)
                    tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
                    tempecModel = setParam(tempecModel,'lb',obj_ec,0.75*gRate_ec);
                    solution_75 = solveLP(tempecModel);
                    if ~isempty(solution_75.x)
                        tempecModel = setParam(ecModel_batch,'obj',ecM_idx,1);
                        tempecModel = setParam(tempecModel,'lb',obj_ec,0.99*gRate_ec);
                        solution_99 = solveLP(tempecModel);
                        if ~isempty(solution_99.x)
                            disp(model.rxnNames{modelIdxs(i)})
                            excRxnIds = [excRxnIds;model.rxnNames(modelIdxs(i))];
                            rates_01  = [rates_01;solution_01.x(ecM_idx)];
                            rates_25  = [rates_25;solution_25.x(ecM_idx)];
                            rates_50  = [rates_50;solution_50.x(ecM_idx)];
                            rates_75  = [rates_75;solution_75.x(ecM_idx)];
                            rates_99  = [rates_99;solution_99.x(ecM_idx)];
                        end
                    end
                end
            end
        end
    end
end
results = table(excRxnIds,rates_01,rates_25,rates_50,rates_75,rates_99); 
sumas = sum(table2array(results(:,2:end)),2);
results = results(sumas>0,:);
writetable(results,'../results/Bmo_production_glc_potential.txt','delimiter','\t','QuoteStrings',false);
    

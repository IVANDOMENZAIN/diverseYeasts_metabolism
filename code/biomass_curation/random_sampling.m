load('../../models/candida_intermedia/cint_GEM_curated.mat')
model = changeMedia_batch(model,'lactose exchange',1);
[solutions, goodRxns] = randomSampling(model,10000);
solutions = full(solutions);
formulas = constructEquations(model);
fluxTable = table(model.rxns,model.rxnNames,model.grRules,formulas,solutions);
writetable(fluxTable,'../../results/randomSampling_WT_lactose.txt','delimiter','\t','QuoteStrings',false)

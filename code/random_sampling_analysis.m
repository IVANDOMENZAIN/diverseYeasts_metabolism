load('../models/candida_intermedia/cint_GEM_curated.mat')
nSamples = 10000;
flux_threshold = 1E-3;
model = setParam(model,'lb','r_4041',0);
model = setParam(model,'ub','r_4041',1000);
model = setParam(model,'obj','r_4041',1);
%run random sampling on wild-type model with unit lactose uptake rate
model = changeMedia_batch(model,'lactose exchange',1);
[fluxTableWT,geneTableWT,geneTable_summaryWT,exch_tableWT] = runRandomSampling(model,nSamples,flux_threshold);
writetable(exch_tableWT,'../results/randomSampling_WT_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
writetable(fluxTableWT,'../results/randomSampling_WT_solutions.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTableWT,'../results/randomSampling_WT_gene_table.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTable_summaryWT,'../results/randomSampling_WT_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)

x = contains(fluxTableWT.Var2,'sorbose exchange');
y = contains(fluxTableWT.Var2,'galactitol exchange');
z = contains(fluxTableWT.Var2,'glucitol exchange');
indexes = [x,y,z];
indexes = sum(indexes,2);
indexes = find(logical(indexes));
newMat = logical(solutions(indexes,:));
indexes = find(sum(newMat,1)>=1);
%indexes = find(sum(logical(),2)>=1);
newT = fluxTableWT(:,1:4);
formulas = constructEquations(model);
newMat = solutions(:,indexes);
indexes = sum(abs(newMat),2)>0;
newMat = newMat(indexes,:);
newTable = table(model.rxns(indexes),model.rxnNames(indexes),formulas(indexes),model.grRules(indexes),newMat);
writetable(newTable,'../results/randomSampling_WT_growth_oxRed_production_solutions.txt','delimiter','\t','QuoteStrings',false)

%set lactose condition
model = changeMedia_batch(model,'lactose exchange',1);
sol = solveLP(model);
model = setParam(model,'lb','r_4041',-0.4999*sol.f);
%model = changeMedia_batch(model,'lactose exchange',1);

[fluxTableWT,geneTableWT,geneTable_summaryWT,exch_tableWT] = runRandomSampling(model,nSamples,flux_threshold);
writetable(exch_tableWT,'../results/randomSampling_WT_growth_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
writetable(fluxTableWT,'../results/randomSampling_WT_growth_solutions.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTableWT,'../results/randomSampling_WT_growth_gene_table.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTable_summaryWT,'../results/randomSampling_WT_growth_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)

%get gal mutant
Galmut = getGALmutant(model);
%set lactose condition
Galmut = changeMedia_batch(Galmut,'lactose exchange',1);
sol = solveLP(Galmut);
Galmut = setParam(Galmut,'lb','r_4041',-0.4999*sol.f);
Galmut = setParam(Galmut,'obj','r_4041',1);
%model = changeMedia_batch(model,'lactose exchange',1);

[fluxTableMut,geneTableMut,geneTable_summaryMut,exch_tableMut] = runRandomSampling(Galmut,nSamples,flux_threshold);
writetable(exch_tableMut,'../results/randomSampling_Mut_growth_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
writetable(fluxTableMut,'../results/randomSampling_Mut_growth_solutions.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTableMut,'../results/randomSampling_Mut_growth_gene_table.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTable_summaryMut,'../results/randomSampling_Mut_growth_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)
%analyse exchanges 
t2 = readtable('../results/randomSampling_Mut_growth_solutions_exch.txt','delimiter','\t');
t1 = readtable('../results/randomSampling_WT_growth_solutions_exch.txt','delimiter','\t');
t1 = table2array(t1(:,5:end));
t2 = table2array(t2(:,5:end));    
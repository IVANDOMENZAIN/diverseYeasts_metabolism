load('../models/candida_intermedia/cint_GEM_curated.mat')
nSamples = 10000;
flux_threshold = 1E-3;
model = setParam(model,'lb','r_4041',0);
model = setParam(model,'ub','r_4041',1000);
model = setParam(model,'obj','r_4041',1);
%run random sampling on wild-type model with unit lactose uptake rate
model = changeMedia_batch(model,'lactose exchange',1);
[fluxTableWT,geneTableWT,geneTable_summaryWT,exch_tableWT] = runRandomSampling(model,nSamples,flux_threshold);
writetable(exch_tableWT,'../results/randomSampling_WT_gal_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
writetable(fluxTableWT,'../results/randomSampling_WT_gal_solutions.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTableWT,'../results/randomSampling_WT_gal_gene_table.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTable_summaryWT,'../results/randomSampling_WT_gal_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)
%set lactose condition
model = changeMedia_batch(model,'lactose exchange',1);
sol = solveLP(model);
model = setParam(model,'lb','r_4041',-0.4999*sol.f);
model = changeMedia_batch(model,'lactose exchange',1);

[fluxTableWT,geneTableWT,geneTable_summaryWT,exch_tableWT] = runRandomSampling(model,nSamples,flux_threshold);
writetable(exch_tableWT,'../results/randomSampling_WT_growth_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
writetable(fluxTableWT,'../results/randomSampling_WT_growth_solutions.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTableWT,'../results/randomSampling_WT_growth_gene_table.txt','delimiter','\t','QuoteStrings',false)
writetable(geneTable_summaryWT,'../results/randomSampling_WT_growth_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)
% 
% %Repeat on GAL mutant
% mutant = getGALmutant(model);
% [fluxTableMut,geneTableMut,geneTable_summaryMut,exch_tableMut] = runRandomSampling(mutant,nSamples,flux_threshold);
% writetable(exch_tableMut,'../../results/randomSampling_GALmut_solutions_exch.txt','delimiter','\t','QuoteStrings',false)
% writetable(fluxTableMut,'../../results/randomSampling_GALmut_solutions.txt','delimiter','\t','QuoteStrings',false)
% writetable(geneTableMut,'../../results/randomSampling_GALmut_gene_table.txt','delimiter','\t','QuoteStrings',false)
% writetable(geneTable_summaryMut,'../../results/randomSampling_GALmut_gene_table_summary.txt','delimiter','\t','QuoteStrings',false)
% 
% 

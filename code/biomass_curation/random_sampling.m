load('../../models/candida_intermedia/cint_GEM_curated.mat')
model = changeMedia_batch(model,'lactose exchange',1);
nSamples = 10000;
[solutions, goodRxns] = randomSampling(model,10000,false,false,true);
solutions = full(solutions);
%get statistical info on flux distributions
mean_flux = mean(solutions,2);
medn_flux = median(solutions,2);
stdd_flux = std(solutions,2);
[grRules,rxnGeneMat,indexes2check] = standardizeGrRules(model);

formulas = constructEquations(model);
sums = sum(abs(solutions),2)>0;
fluxTable = table(model.rxns(sums),model.rxnNames(sums),model.grRules(sums),formulas(sums),solutions(sums,:));
writetable(fluxTable,'../../results/randomSampling_WT_lactose.txt','delimiter','\t','QuoteStrings',false)
%Get a logical variant of the solutions matrix
sol_logical = logical(solutions);
sol_active = sum(sol_logical,2);
%visualize the number of times that reactions are active 
histogram(sol_active)
 xlabel('reaction occurrences')
 ylabel('frequency')
flux_threshold = 1E-2;
sol_logical_stringent = abs(solutions)>=flux_threshold;
sol_active_stringent = sum(sol_logical_stringent,2);
%visualize the number of times that reactions are active (stringent)
histogram(sol_active_stringent)
xlabel('reaction occurrences')
ylabel('frequency')
%focus on highly active reactions
occurrence_threshold = 0;
high_occurr_rxns = find(sol_active_stringent>occurrence_threshold*nSamples);
high_occurr_table = solutions(high_occurr_rxns,:);
mean_flux = mean(solutions,2);
medn_flux = median(solutions,2);
stdd_flux = std(solutions,0,2);
mode_flux = mode(solutions,2);
reduced_fluxTable = table(model.rxns(high_occurr_rxns),model.rxnNames(high_occurr_rxns),model.grRules(high_occurr_rxns),formulas(high_occurr_rxns),sol_active_stringent(high_occurr_rxns),mean_flux(high_occurr_rxns),medn_flux(high_occurr_rxns),stdd_flux(high_occurr_rxns),mode_flux(high_occurr_rxns));
%let's get a gene-solutions matrix
G= numel(model.genes);
S = nSamples;
R = numel(model.rxns);
geneRndMat = zeros(G,S);
model.rxnGeneMat = rxnGeneMat;
newMat = rxnGeneMat';
temp = abs(solutions)>=flux_threshold;
newMat = logical(newMat*(temp));
%
sum_genes = sum(newMat,2);
mean_genes = mean(newMat,2);
medn_genes= median(newMat,2);
stdd_genes = std(newMat,0,2);
mode_genes = mode(newMat,2);
geneTable = table(model.genes,model.geneShortNames,model.proteins,sum_genes,mean_genes,medn_genes,stdd_genes,mode_genes);
%isolate the lactose-related genes
lac_genes = {'Seq_2552' 'Seq_2479' 'Seq_3332' 'xyl1_2' 'xyl1' 'xyl1_3' 'Seq_2272' 'Seq_2189' 'Seq_4294' 'Seq_1935' 'Seq_3460' 'Seq_5357' 'Seq_2923'};
[a,b] = ismember(lac_genes,model.genes);
lac_gene_table = geneTable(b',:);
oxRed_genes = {'Seq_2552' 'xyl1_2' 'xyl1' 'xyl1_3' 'Seq_2272' 'Seq_2189' 'Seq_4294' 'Seq_5357' 'Seq_2923'};
[a,b] = ismember(oxRed_genes,model.genes);
oxred_gene_table = geneTable(b',:);
oxred_gene_solutions = newMat(b,:);
sum_oxred_simul = sum(oxred_gene_solutions,1);
all_oxred = find(sum_oxred_simul>=2);
reduced_solutions_AlloxredGenes = solutions(:,all_oxred); 
allOxRed_solutions = table(model.rxns,model.rxnNames,model.grRules,formulas,reduced_solutions_AlloxredGenes);
writetable(allOxRed_solutions,'../../results/randomSampling_WT_lactose_oxredGenes.txt','delimiter','\t','QuoteStrings',false)
[a,b] = getExchangeRxns(model);
allOxRed_solutions_exch = allOxRed_solutions(b,:);
%identify secreted products
logicalSol = (reduced_solutions_AlloxredGenes(b,:)>0);
sumLogicalSol = sum(logicalSol,2);
[a,b] = sort(sumLogicalSol,'descend');
allOxRed_solutions_exch = allOxRed_solutions_exch(b,:);
%get a table that sorts all solutions by number of times metabolites are secreted
topMets = find(a>=5);
allOxRed_solutions_exch = allOxRed_solutions_exch(topMets,:);
%interesting mets 
target_mets = [allOxRed_solutions_exch.Var1,allOxRed_solutions_exch.Var4];
allOxRed_solutions_exch = allOxRed_solutions_exch(5:end,:);
%now we got a list of the interesting exchange reactions to explore in the
%solutions from random sampling
%let's try to understand how is L-sorbose secreted by the model
x=find(strcmp(model.rxns,'r_1909'));
indexes = find(reduced_solutions_AlloxredGenes(x,:)>flux_threshold);
sorbose_solutions = table(model.rxns,model.rxnNames,model.grRules,formulas,reduced_solutions_AlloxredGenes(:,indexes));


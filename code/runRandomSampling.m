function [fluxTable,geneTable,geneTable_summary,exch_table] = runRandomSampling(model,nSamples,flux_threshold)
solutions  =randomSampling(model,nSamples,true,false,true,[],true);
solutions = full(solutions);
%get statistical info on flux distributions
formulas  = constructEquations(model);
fluxTable = table(model.rxns,model.rxnNames,model.grRules,formulas,solutions);
%Get a logical variant of the solutions matrix
sol_logical = logical(solutions);
sol_active  = sum(sol_logical,2);
%let's get a gene-solutions matrix
[~,rxnGeneMat] = standardizeGrRules(model);
G = numel(model.genes);
R = numel(model.rxns);
geneRndMat       = zeros(G,nSamples);
model.rxnGeneMat = rxnGeneMat;
newMat = rxnGeneMat';
temp   = abs(solutions) >= flux_threshold;
newMat = logical(newMat*temp);
%
sum_genes  = sum(newMat,2);
mean_genes = mean(newMat,2);
medn_genes = median(newMat,2);
stdd_genes = std(newMat,0,2);
mode_genes = mode(newMat,2);
if length(model.proteins)== length(model.genes)
    geneTable_summary  = table(model.genes,model.geneShortNames,model.proteins,sum_genes,mean_genes,medn_genes,stdd_genes,mode_genes);
    geneTable  = table(model.genes,model.geneShortNames,model.proteins,newMat);
else
    geneTable_summary  = table(model.genes,model.geneShortNames,sum_genes,mean_genes,medn_genes,stdd_genes,mode_genes);
    geneTable  = table(model.genes,model.geneShortNames,newMat);    
end
%get a secretions table
[~,b] = getExchangeRxns(model);
solutions_exch = solutions(b,:);
exch_table = table(model.rxns(b),model.rxnNames(b),model.grRules(b),formulas(b),solutions_exch);
%identify secreted products
logicalSol    = (solutions_exch>0);
sumLogicalSol = sum(logicalSol,2);
%get a table that sorts all solutions by number of times metabolites are secreted
[~,b] = sort(sumLogicalSol,'descend');
exch_table = exch_table(b,:);
end
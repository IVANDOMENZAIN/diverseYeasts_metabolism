function [coeff,score,latent,tsquared,explained]= PCA_analysis(trackedElement,model,matrix,field)
if strcmp(field,'met')
    x = find(strcmp(model.metNames,trackedElement));
    x = x(1);
    [matrix] = getMetTurnOver(matrix,model);
elseif strcmp(field,'rxn')
    x = find(strcmp(model.rxns,trackedElement));
end
vector = matrix(x,:)/max(matrix(x,:));
colors = zeros(2*nSamples,3);
colors(:,1) = vector;
colors(:,2) = 0.4;
colors(:,3) = 0.5;
[coeff,score,latent,tsquared,explained] = pca(matrix');
PC1 = score(:,1);
PC2 = score(:,2);
x_lab = ['PC1: ' num2str(round(explained(1),2)) '% of variance'];
y_lab = ['PC2: ' num2str(round(explained(2),2)) '% of variance'];
fig = scatter(PC1,PC2,40,colors,'fill');
xlabel(x_lab,"FontSize",18);
ylabel(y_lab,"FontSize",18);
title(['PCA:  ' trackedElement  ' turnover rate'],"FontSize",20)
saveas(fig,['../results/figures/' trackedElement '_' type '_PCA_wtMut_growth_lactose.fig'])
saveas(fig,['../results/figures/' trackedElement '_' type '_PCA_wtMut_growth_lactose.pdf'])
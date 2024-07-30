% function [coeff,score,latent,tsquared,explained]= PCA_analysis(trackedElement,model,matrix,field)
bckgnd = 'WT';
csoure = 'gal';
condition = 'clim';
nSamples = 10000;

t = readtable(['../results/randomSampling_' bckgnd '_' csoure '_growth_' condition '_solutions.txt'],'delimiter','\t');
t = t(:,5:end);
t = table2array(t);
load('../models/candida_intermedia/cintGEM_curated.mat')
[matrix] = getMetTurnOver(t,model);
%matrix = t;
[coeff,score,latent,tsquared,explained] = pca(matrix');

close all
hold off
set(gca,'FontSize',20, 'FontName', 'Courier')
trackedElement = 'ATP';
x = find(strcmp(model.metNames,trackedElement));
x = x(model.metComps(x)==1);
vector = matrix(x,:)/(max(matrix(x,:)));
colors = zeros(nSamples,3);
colors(:,1) = vector;
colors(:,2) = 0.8*vector;
colors(:,3) = 1-vector;
PC1 = score(:,1);
PC2 = score(:,2);

x_lab = ['PC1: ' num2str(round(explained(1),2)) '% of variance'];
y_lab = ['PC2: ' num2str(round(explained(2),2)) '% of variance'];
fig = scatter(PC1,PC2,200,colors,'LineWidth',3.5,'Marker','o');
%title(['PCA:  ' trackedElement  ' turnover rate'],"FontSize",20)
xlabel(x_lab,"FontSize",20);
ylabel(y_lab,"FontSize",20);
%set(gca,'FontSize',18)
hold on

x = find(strcmp(model.metNames,'alpha-D-galactose 1-phosphate'));
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>0);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0;
colors(:,2) = 0;%vector(idx);
colors(:,3) = 0;%vector(idx);
fig = scatter(PC1(idx),PC2(idx),90,colors,'LineWidth',2.5,'Marker','diamond');
hold on

x = find(strcmp(model.metNames,'galactitol'));
x = x(1);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>0);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0;
colors(:,2) = 1;%vector(idx);
colors(:,3) = 1;%vector(idx);
fig = scatter(PC1(idx),PC2(idx),80,colors,'LineWidth',1,'Marker','square');
hold on

% x = find(strcmp(model.rxnNames,'aldose reductase (NAPDH)'));
% vector = t(x,:)/max(t(x,:));
% idx = find(vector>0);
% colors = zeros(numel(idx),3)+1;
% colors(:,1) = 0;
% colors(:,2) = 1;%vector(idx);
% colors(:,3) = 1;%vector(idx);
% fig = scatter(PC1(idx),PC2(idx),80,colors,'filled','LineWidth',1.5,'Marker','square');
% hold on

x = find(strcmp(model.metNames,'biomass'));
x = x(1);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.9);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 1;%vector(idx);
colors(:,2) = 0;
colors(:,3) = 0.15;
alpha = vector(idx);
fig = scatter(PC1(idx),PC2(idx),150,colors,'filled','LineWidth',1.5,'Marker','v');
hold on
set(gca,'FontSize',18)

%xlim([-300 300])
% %ylim([-100 200])
% %saveas(fig,['../results/figures/' trackedElement '_met_PCA_' bckgnd '_growth_' csoure '_' condition '.fig'])
% %saveas(fig,['../results/figures/' trackedElement '_met_PCA_' bckgnd '_growth_' csoure '_' condition '.pdf'])
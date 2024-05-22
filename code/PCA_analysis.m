% function [coeff,score,latent,tsquared,explained]= PCA_analysis(trackedElement,model,matrix,field)
bckgnd = 'WT';
csoure = 'gal';
condition = 'oxlim';
nSamples = 10000;


t = readtable(['../results/randomSampling_' bckgnd '_' csoure '_growth_' condition '_solutions.txt'],'delimiter','\t');
t = t(:,5:end);
t = table2array(t);
load('../models/candida_intermedia/cintGEM_curated.mat')
[matrix] = getMetTurnOver(t,model);
[coeff,score,latent,tsquared,explained] = pca(matrix');

hold off
set(gca,'FontSize',20, 'FontName', 'Courier')
trackedElement = 'ATP';
x = find(strcmp(model.metNames,trackedElement));
x = x(1);
vector = sum(abs(t),1);
vector = vector/max(vector);
vector = matrix(x,:)/(max(matrix(x,:)));
%vector = vector';
%vector(vector>1)=1;
colors = zeros(nSamples,3);
% colors(1:10000,3) = 1;
% colors(:,1) = 0.9*vector;
% colors(:,2) = 0.75*vector;

colors(:,1) = vector;
colors(:,2) = 0.9*vector;
colors(:,3) = 1-vector;
PC1 = score(:,1);
PC2 = score(:,2);

x_lab = ['PC1: ' num2str(round(explained(1),2)) '% of variance'];
y_lab = ['PC2: ' num2str(round(explained(2),2)) '% of variance'];
fig = scatter(PC1,PC2,100,colors,'LineWidth',3);
%title(['PCA:  ' trackedElement  ' turnover rate'],"FontSize",20)

hold on

x = find(strcmp(model.metNames,'galactitol'));
x = x(2);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.5);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0;
colors(:,2) = 1;%vector(idx);
colors(:,3) = 1;%vector(idx);
fig = scatter(PC1(idx),PC2(idx),40,colors,'filled','LineWidth',1,'Marker','square');
% 
hold on




hold on
x = find(strcmp(model.metNames,'D-tagatose'));
x = x(2);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.1);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 1;
colors(:,2) = 0.4;
colors(:,3) = 1;
fig = scatter(PC1(idx),PC2(idx),200,colors,'filled','LineWidth',1,'Marker','pentagram');



x = find(strcmp(model.metNames,'L-sorbose'));
x = x(2);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.1);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0.6;
colors(:,2) = 0.6;
colors(:,3) = 0.6;
fig = scatter(PC1(idx),PC2(idx),100,colors,'filled','LineWidth',1,'Marker','v');

hold on
x = find(strcmp(model.metNames,'D-glucitol'));
x = x(2);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.1);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0;
colors(:,2) = 0;
colors(:,3) = 0;
fig = scatter(PC1(idx),PC2(idx),100,colors,'filled','LineWidth',1,'Marker','o');

hold on
x = find(strcmp(model.metNames,'biomass'));
x = x(1);
vector = matrix(x,:)/max(matrix(x,:));
idx = find(vector>=0.99);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 1;%vector(idx);
colors(:,2) = 0;
colors(:,3) = 0.15;
alpha = vector(idx);
fig = scatter(PC1(idx),PC2(idx),100,colors,'LineWidth',1.5,'Marker','*');
set(gca,'FontSize',18)
xlabel(x_lab,"FontSize",20);
ylabel(y_lab,"FontSize",20);
%xlim([-300 300])
%ylim([-100 200])
saveas(fig,['../results/figures/' trackedElement '_met_PCA_' bckgnd '_growth_' csoure '_' condition '.fig'])
saveas(fig,['../results/figures/' trackedElement '_met_PCA_' bckgnd '_growth_' csoure '_' condition '.pdf'])

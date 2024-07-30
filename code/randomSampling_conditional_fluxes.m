csources = {'lac' 'gal'};
conditions = {'clim' 'oxlim'};
nSamples  = 10000;

load('../models/candida_intermedia/cintGEM_curated.mat')
ethanol_exh    = find(strcmp(model.rxnNames,'ethanol exchange'));
galactitol_exh = find(strcmp(model.rxnNames,'galactitol exchange'));
growth_idx     = find(strcmp(model.rxnNames,'growth'));
oxy_exh        = find(strcmp(model.rxnNames,'oxygen exchange'));

galIdxs = {};
fluxDists = [];
fluxDists_galactitol =[];
fluxDists_ethanol =[];
for csource = {'lac' 'gal'}
    for condition = {'clim' 'oxlim'}
        condi = condition{1};
        source = csource{1};
        t = readtable(['../results/randomSampling_WT_' source '_growth_' condi '_solutions.txt'],'delimiter','\t');
        annotatio = t(:,1:4);
        t = t(:,5:end);
        t = table2array(t);
        %
        gidx = find(t(galactitol_exh,:)>0);
        eidx = find(t(ethanol_exh,:)>0);
        just_galactitol = setdiff(gidx,eidx);
        just_ethanol = setdiff(eidx,gidx);
        if strcmp(condi,'clim')
            galIdxs = [galIdxs,{just_galactitol}];
            fluxDists_galactitol = [fluxDists_galactitol,t(:,just_galactitol)];
            fluxDists_ethanol    = [fluxDists_ethanol,t(:,just_ethanol)];
        end
        disp([condi ' / ' source ' / There are ' num2str(numel(just_galactitol)) ' flux distributions producing galactitol but not ethanol'])
        disp([condi ' / ' source ' / There are ' num2str(numel(just_ethanol)) ' flux distributions producing ethanol but not galactitol'])

    end
end


t1 = array2table(fluxDists_galactitol);
t2 = array2table(fluxDists_ethanol);
t3 = table(annotatio,t1,t2);

newMat = ([fluxDists_galactitol,fluxDists_ethanol]);
[~,nSamples] = size(newMat);

[coeff,score,latent,tsquared,explained] = pca(newMat');



uno = 5;
dos = 6; 
PC1 = score(:,uno);
PC2 = score(:,dos);

hold off
x = find(strcmp(model.rxnNames,'growth'));
%x = x(model.metComps(x)==1);
vector = newMat(x,:)/(max(newMat(x,:)));
colors = zeros(nSamples,3);
colors(:,1) = vector;
colors(:,2) = 0.8*vector;
colors(:,3) = 1-vector;
x_lab = ['PC1: ' num2str(round(explained(uno),2)) '% of variance'];
y_lab = ['PC2: ' num2str(round(explained(dos),2)) '% of variance'];
fig = scatter(PC1,PC2,200,colors,'LineWidth',3.5,'Marker','o');
%title(['PCA:  ' trackedElement  ' turnover rate'],"FontSize",20)
xlabel(x_lab,"FontSize",20);
ylabel(y_lab,"FontSize",20);
hold on

x = find(strcmp(model.rxnNames,'ethanol exchange'));
vector = newMat(x,:)/max(newMat(x,:));
idx = find(vector>0);

%idx = find(vector>=0.5);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 0;
colors(:,2) = 0;
colors(:,3) = 0;
fig = scatter(PC1(idx),PC2(idx),300,colors,'filled','LineWidth',1,'Marker','pentagram');
hold on

x = find(strcmp(model.rxnNames,'galactitol exchange'));
vector = newMat(x,:)/max(newMat(x,:));
idx = find(vector>0);
%idx = find(vector>=0.5);
colors = zeros(numel(idx),3)+1;
colors(:,1) = 1;
colors(:,2) = 0;
colors(:,3) = 0;
fig = scatter(PC1(idx),PC2(idx),300,colors,'filled','LineWidth',1,'Marker','pentagram');
hold on








bckgnds = {'Mut'};
csoure = 'lac';
%conditions = {'clim' 'oxlim' 'nlim'};
conditions = {'clim'};
nSamples = 10000;
flux_threshold = 0.01;

load('../models/candida_intermedia/cintGEM_curated.mat')

indexes(1) =  find(strcmp(model.rxnNames,'galactitol exchange'));
indexes(2) =  find(strcmp(model.rxnNames,'L-sorbose exchange'));
indexes(3) =  find(strcmp(model.rxnNames,'D-glucitol exchange'));
indexes(4) =  find(strcmp(model.rxnNames,'D-tagatose exchange'));
indexes(5) =  find(strcmp(model.rxnNames,'ethanol exchange'));
[a,indexes] = getExchangeRxns(model);

secretions = zeros(numel(bckgnds)*numel(conditions),numel(indexes));
averages = secretions;
counter = 1;
tags = {};
for i=1:numel(bckgnds)
    bckgnd = bckgnds{i};
    for j=1:numel(conditions)
        condition = conditions{j};
        t = readtable(['../results/randomSampling_' bckgnd '_' csoure '_growth_' condition '_solutions.txt'],'delimiter','\t');
        t = t(:,5:end);
        t = table2array(t);       
        for k=1:numel(indexes)
            positions = find(t(indexes(k),:)>flux_threshold);
            secretions(counter,k) = sum(t(indexes(k),:)>flux_threshold);
            averages(counter,k) = mean(t(indexes(k),positions));
        end
        tags{counter} = [bckgnd '_' condition];
        counter = counter + 1;
    end
end
rxnNames   = model.rxnNames(indexes);
ocurrences = sum(secretions,1);
[a,b] = sort(ocurrences,'descend');

secretions = array2table(secretions);
averages   = array2table(averages);

secretions.Properties.VariableNames = rxnNames;
secretions.Properties.RowNames = tags;

averages.Properties.VariableNames = rxnNames;
averages.Properties.RowNames = tags;

secretions = rows2vars(secretions);
averages = rows2vars(averages);

secretions = secretions(b,:);
averages = averages(b,:);


writetable(secretions,['../results/randomSampling_Mut_' csoure '_growth_secretions.txt'],'delimiter','\t','QuoteStrings',false)
writetable(averages,['../results/randomSampling_Mut_' csoure '_growth_Average_secretion_flux.txt'],'delimiter','\t','QuoteStrings',false)


%%% Script adapted from KMGO60 Systems Biology course
%  RNAseq DE analysis
%
%  Christoph Boerlin
%  Ivan Domenzain
%
%  Last edited. 2019-09-17

%%%Step 1: Load RNAseq Data and normalize the data

%Load RNAseq data
rawCounts = readtable('../../data/RNAseq/normalized_counts.txt','delimiter','\t');

%Create logical vectors for easy access of the different conditions
%now rawCounts(:,HiT) outputs only the three columns with the data for
%the high temperature _Csrc condition.
samples = rawCounts(:,3:end).Properties.VariableNames;
glu = startsWith(samples,'G20_');
cel = startsWith(samples,'C20_');
gal = startsWith(samples,'GA20_');
lac = startsWith(samples,'L20_');
xyl = startsWith(samples,'X20_');

samplesPos = [glu; cel; gal; lac; xyl];

conditions = {'glu' 'cel' 'gal' 'lac' 'xyl'};
%Get only the counts from the RNAseq data
rawCountsNum = rawCounts{:,3:end};
%Let's plot all the reads in the samples as boxplots
boxplot(rawCountsNum)
title('Raw counts')
ylabel('Read counts')
xlabel('Samples')
%Samples are pretty alligned, the counts values are remarkably low, this
%suggests a possible prior log2 transformation, let's recover the original
%values
rawCountsLog = rawCountsNum;
rawCountsNum = (2.^rawCountsNum);
%Estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(rawCountsNum,2);
%Data normalization
%First, get the positions for all of the genes with non-zero reads
nz = pseudoRefSample > 0;
%Get ratios of expression, dividing each read by its corresponding gene
%average expression across samples
ratios = rawCountsNum(nz,:)./pseudoRefSample(nz);
%Get sizeFactors for normalization of each column (library of gene reads) 
sizeFactors = median(ratios,1);
% normalize the raw counts using the calculated normalization factors
normCounts = rawCountsNum./sizeFactors;
%Let's take a look to the normalization effects on the gene counts
%distributions. Generate boxplots for each of the samples distributions in
%both rawCountsNum and also normCounts
subplot(1,2,1);
boxplot(rawCountsNum)
title('Raw counts')
ylabel('Read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(normCounts)
title('Normalized counts')
ylabel('Read counts')
xlabel('Samples')

% Lets visualize this on log2 scale
subplot(1,2,1);
boxplot(log2(rawCountsNum))
title('Raw counts')
ylabel('log_{2} read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(log2(normCounts))
title('Normalized counts')
ylabel('log_{2} read counts')
xlabel('Samples')

%%%Step 2: PCA analysis

%Perform PCA
[pc, zscores, pcvars, ~, explained] = pca(normCounts');

%Pick 5 colors to distinguish the 5 different conditions. Colormap returns
%the a mapped normalized RGB values for the specified amount of colorlevels
%in a given palette (jet in this case)
colorScheme = colormap(jet(5));
hold on
%Get a PCA plot showing PC1 and PC2 for all the samples in the 5 conditions
% **Note: PC's are stored in the zscores table in which for each sample 
%         (columns) the 15 first principal components are shown (rows).
condStr = {};
for i=1:length(conditions)
    %Get the condition name
    nameTMP = conditions{i};
    PC1     = zscores(find(samplesPos(i,:)),1);
    PC2     = zscores(find(samplesPos(i,:)),2);
    condStr = [condStr;nameTMP];
    scatter(PC1,PC2,50,colorScheme(i,:),'fill');
end
xlabel(['First Principal Component ' num2str(explained(1)) '%']);
ylabel(['Second Principal Component ' num2str(explained(2)) '%']);
title('Principal Component Scatter Plot');
legend(condStr)
hold off

%%%Step 3: Test for differential gene expression

%pick one of the four _Csrc condition that you would like to compare to
%the control
for j = 1:length(conditions)
    cond = conditions{j};
    if ~strcmpi(cond,'glu')
        disp(cond)
        normCounts_Csrc  = normCounts(:,samplesPos(j,:));
        normCountsref    = normCounts(:,samplesPos(1,:));
        
        %Perform the test for differential expression using a negative binomial model
        tLocal = nbintest(normCountsref,normCounts_Csrc,'VarianceLink','LocalRegression');
        % correct for multiple testing using the Benjamini Hochberg correction
        padj = mafdr(tLocal.pValue,'BHFDR',true);
        
        %create overview table of the results including the mean for both
        %conditions, the resulting log2 fold change and the calculated pValue
        meanRef    = mean(normCountsref,2);
        mean_Csrc = mean(normCounts_Csrc,2);
        log2FC     = log2(mean_Csrc./meanRef);
        geneTable  = table(rawCounts.geneNames, meanRef,mean_Csrc,log2FC,tLocal.pValue,padj);
        %Add row and column names
        geneTable.Properties.RowNames      = rawCounts.genes;
        geneTable.Properties.VariableNames = {'geneNames' 'Mean_Ref','Mean__Csrc','Log2_FC','pVal','adjPVal'};
        
        %create a Volcano Plot for visual inspection of the data
        thresholds = [2 3];
        
        scatter(geneTable.Log2_FC,-log10(geneTable.adjPVal),30,'fill')
        DEgenesPos = abs(geneTable.Log2_FC)>thresholds(1) & -log10(geneTable.adjPVal)>thresholds(2);
        DEgenesDown = (geneTable.Log2_FC)<-thresholds(1) & -log10(geneTable.adjPVal)>thresholds(2);
        DEgenesUp = (geneTable.Log2_FC)>thresholds(1) & -log10(geneTable.adjPVal)>thresholds(2);
        disp(['There are ' num2str(sum(DEgenesPos)) ' DE genes, from which ' num2str(sum(DEgenesDown)) ' are dReg and ' num2str(sum(DEgenesUp)) ' are upReg'])
        disp(' ')
        scatter(geneTable.Log2_FC,-log10(geneTable.adjPVal),30,'fill')
        mkdir('../results/RNA_DE_analysis')
        writetable(geneTable,['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' cond '.txt'],'delimiter','\t','QuoteStrings',false,'WriteRowNames',true)
    end
end
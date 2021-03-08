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
% 
% %%%Step 4: First analysis of DE genes
% 
% %Manual inspection of genes with highest differential expression
% geneTable = sortrows(geneTable,'adjPVal','ascend');
% disp(geneTable(1:10,:))
% 
% %##################
% %TASK: 
% %Do you find any interesting genes or patterns in this top differentially expressed genes?
% %Combine the geneTable with the geneDescription table (from Data/GeneDescriptions.csv) for getting a short
% %description for each gene.
% %
% %HINT:
% %Use the intersection function to find the rows corresponding to the same gene
% %##################
% %%%Step 5: Find Associated GO Terms
% 
% %load assocition of GO Terms and genes
% GoTerms_table = readtable('../data/GoTermsMapping.txt','delimiter','\t');
% GoTermsIDs    = unique(GoTerms_table.GoTerm);
% 
% %convert it to a map (the matlab variant of a dictionary) for easier access
% %A map consists of unique key - value pairs and has the advantage that with
% %the key you can easily access the corresponding value. If you have a
% %map called "Universities" and you want to look up the value associated to 
% %the Key "Chalmers" use Universities('Chalmers') and you would get the 
% %data associated to Chalmers.
% GoTermsGeneMap = containers.Map();
% for i = 1:height(GoTerms_table)
%     key = GoTerms_table{i,'GeneName'}{1};
%     if isKey(GoTermsGeneMap,key)
%         GoTermsGeneMap(key) = [GoTermsGeneMap(key),GoTerms_table{i,'GoTerm'}{1}];
%     else
%         GoTermsGeneMap(key) = {GoTerms_table{i,'GoTerm'}{1}};
%     end
% end
% fprintf('Number of annotated genes related to functional process is %d.\n',GoTermsGeneMap.Count)
% fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique(GoTerms_table.GoTerm)))
% fprintf('Number of gene-GO term associations is %d.\n',numel(GoTerms_table))
% 
% %Look up associated GO Terms for the top differentially expressed gene
% selectedGene      = geneTable.Properties.RowNames{1};
% associatedGoTerms = GoTermsGeneMap(selectedGene);
% disp(['Associated GO Terms for gene ',selectedGene])
% disp(associatedGoTerms)
% 
% %The GO Term IDs are not informative, therefore we need to load the
% %descriptions for them
% GO = geneont('File','../data/GoTerms.obo');
% 
% %Look up the details of the first GO Term
% %The lookup function needs just the ID number (in number format) as input
% GoTermID         = str2double(associatedGoTerms{1}(4:end));
% GoTermName       = GO(GoTermID).Terms.Name;
% GoTermDefinition = GO(GoTermID).Terms.Definition;
% %output result
% disp(['GO Term ',num2str(GoTermID),' - ',GoTermName,' : ',GoTermDefinition])
% 
% %##################
% % TASK: 
% % Create a loop to print the details for each of the associated GO Terms
% % for the top1 differentially expressed gene.
% %
% % What can you learn from it? Compare the knowledge you got with the
% % summary paragraph on the Saccharomyces Genome Database (yeastgenome.org)
% %##################
% 
% %%%Step 6: Find enriched GO Terms in differentially expressed genes
% 
% %split gene list into DE and non DE genes
% indexDE       = geneTable.adjPVal<=0.01 & abs(geneTable.Log2_FC)>=2;
% geneListDE    = geneTable.Properties.RowNames(indexDE);
% geneListNonDE = geneTable.Properties.RowNames(~indexDE);
% %Select GO Term, in this case the top1 is chosen
% rankedPosition = 1;
% GoTermID       = str2double(associatedGoTerms{rankedPosition}(4:end));
% %GO terms full IDs are of the form GO:XXXXXXX
% full_GoTerm_Id = sprintf('GO:%07d',GoTermID);
% 
% %################## 'Alternative #1 
% %Check enrichment for selected GO Terms
% %Count number of occurences for the GO Term in all genes and in all DE
% %genes
% GoTermCountAll = 0;
% GoTermCountDE  = 0;
% for key = GoTermsGeneMap.keys
%     associatedGoTerms = (GoTermsGeneMap(key{1}));
%     %Search the GO term ID in the associated GO terms cell array
%     if any(strcmp(associatedGoTerms,full_GoTerm_Id))
%         %
%         if any(strcmp(geneListDE,key{1}))
%             GoTermCountDE = GoTermCountDE+1;
%         end
%         GoTermCountAll = GoTermCountAll+1;
%     end
% end
% %################## Alternative 2
% GoTermsIDs_table  = GoTerms_table.GoTerm;
% GoTermCountAll = 0;
% GoTermCountDE  = 0;
% presence       = find(strcmp(GoTermsIDs_table,full_GoTerm_Id));
% if ~isempty(presence)
%     relatedGenes  = GoTerms_table.GeneName(presence);
%     [~,indexesDE] = intersect(geneListDE,relatedGenes);
%     GoTermCountDE  = length(indexesDE);
%     GoTermCountAll = length(presence);
% end
% %##################
% 
% disp(['GO Term ',num2str(GoTermID),' is in ',num2str(GoTermCountDE),...
%     ' out of ',num2str(numel(geneListDE)),' DE Genes and in total there are ',...
%     num2str(GoTermCountAll),' occurences in all ',num2str(numel(GoTermsGeneMap.keys)),' Genes'])
% 
% %run hypergeometric test to assess significance
% pHyperGeo=hygepdf(GoTermCountDE,numel(GoTermsGeneMap.keys),GoTermCountAll,numel(geneListDE));
% disp(['The probablity for this is ',num2str(pHyperGeo)])
% %##################
% % TASK: 
% % Create a loop that goes over all existing GoTerms (GoTermsIDs) and
% % calculates the enrichment. Don't forget to use a correction for multiple
% % testing, e.g. the Benjamini Hochberg Correction as for the DE detection
% % (mafdr).
% % Sort GO Terms by adjusted pValue and print out every one that has a
% % pValue <=0.01
% %##################
% 
% 

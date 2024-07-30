%function [pearson_mat,pValues_mat] = gene_exp_corr(dataTable)
dataTable = readtable('../data/RNAseq/normalized_counts.txt','Delimiter','\t');
countMatrix = table2array(dataTable(:,3:end));
countMatrix = 2.^countMatrix;

%[a,b] = ismember(lacGenes,dataTable.genes);
%countMatrix = countMatrix(b,:);

minVals = min(countMatrix,[],1);
[nGenes,nSamples] = size(countMatrix);
pearson_mat = zeros(nGenes,nGenes);
pValues_mat = zeros(nGenes,nGenes);
for i=1:nGenes
    vector1 = countMatrix(i,:);
    parfor j=1:nGenes
        if j==i
            pearson_mat(i,j) = 1;
        elseif j<i
            vector2 = countMatrix(j,:);
            [cmat,pVal] = corrcoef(vector1,vector2);
            pearson_mat(i,j) = cmat(2,1);
            pValues_mat(i,j) = -log10(pVal(2,1));
        end
    end
end

for i=1:nGenes
    for j=1:nGenes
        if j>i
            pearson_mat(i,j) = pearson_mat(j,i);
            pValues_mat(i,j) = pValues_mat(j,i);
        end
    end
end
%end
lacGenes = {'CIC11T00000000334' 'CIC11T00000000893' 'CIC11T00000005922' 'CIC11T00000003459' ...
            'CIC11T00000002481' 'CIC11T00000001944' 'CIC11T00000001159' 'CIC11T00000000932' ...
            'CIC11T00000000159' 'CIC11T00000003249' 'CIC11T00000000287' 'CIC11T00000005567' ...
            'CIC11T00000002104' 'CIC11T00000005750' 'CIC11T00000001388' 'CIC11T00000003915' 'CIC11T00000004651'};
lacIDs   = {'xyl1' 'xyl1-2' 'xyl1-3' 'lad' ...
            'lxr' 'xyl2' 'Seq1183' 'Seq4936' ...
            'gal10' 'gal10-2' 'gal1' 'gal1-2' ...
            'gal7' 'LAC4' 'LAC12_3' 'LAC9' 'LAC9_2'};

lacIDs = upper(lacIDs);
lacIDs = strrep(lacIDs,'_','-')
[a,b] = ismember(lacGenes,dataTable.genes);
lacMatrixRo = pearson_mat(b,b);
newVect = reshape(pearson_mat,nGenes*nGenes,1);

figure
lacMatrixRo = round(lacMatrixRo,2);
h = heatmap(lacMatrixRo,'MissingDataColor','w')
h.FontSize = 18
h.XDisplayLabels = lacIDs
h.YDisplayLabels = lacIDs

figure
lacMatrixPV = pValues_mat(b,b);
lacMatrixPV = round(lacMatrixPV,2);
lacMatrixPV(lacMatrixPV>3) =3;
h = heatmap(lacMatrixPV,'MissingDataColor','w')
h.FontSize = 18
h.XDisplayLabels = lacIDs
h.YDisplayLabels = lacIDs
hold on
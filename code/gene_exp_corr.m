dataTable = readtable('../data/RNAseq/normalized_counts.txt','Delimiter','\t');
countMatrix = table2array(dataTable(:,3:end));
countMatrix = 2.^countMatrix;

lacGenes = {'xyl1' 'xyl1_2' 'xyl1_3' 'Seq_2552' 'Seq_2272' 'Seq_2189' 'Seq_2479' 'Seq_3332' 'Seq_1935' 'Seq_4294' 'Seq_3460'};
lacGenes = {'CIC11T00000000334' 'CIC11T00000000893' 'CIC11T00000005922' 'CIC11T00000003459' 'CIC11T00000002481' 'CIC11T00000001944' 'CIC11T00000000159' 'CIC11T00000003249' 'CIC11T00000000287' 'CIC11T00000005567' 'CIC11T00000002104'};
%[a,b] = ismember(lacGenes,dataTable.genes);
%countMatrix = countMatrix(b,:);

minVals = min(countMatrix,[],1);
[nGenes,nSamples] = size(countMatrix);
pearson_mat = zeros(nGenes,nGenes);
for i=1:nGenes
    for j=1:nGenes
        if j==i
           pearson_mat(i,j) = 1;
        elseif j<i
            vector1 = countMatrix(i,:);
            vector2 = countMatrix(j,:);
            cmat = corrcoef(vector1,vector2);
            pearson_mat(i,j) = cmat(2,1);
        end
    end
end

for i=1:nGenes
    for j=1:nGenes
        if j>i
            pearson_mat(i,j) = pearson_mat(j,i);
        end
    end
end
lacGenes = {'xyl1' 'xyl1_2' 'xyl1_3' 'Seq_2552' 'Seq_2272' 'Seq_2189' 'Seq_2479' 'Seq_3332' 'Seq_1935' 'Seq_4294' 'Seq_3460'};
lacGenes = {'CIC11T00000000334' 'CIC11T00000000893' 'CIC11T00000005922' 'CIC11T00000003459' 'CIC11T00000002481' 'CIC11T00000001944' 'CIC11T00000001159' 'CIC11T00000000932' 'CIC11T00000000159' 'CIC11T00000003249' 'CIC11T00000000287' 'CIC11T00000005567' 'CIC11T00000002104'};
lacIDs   = {'xyl1' 'xyl1-2' 'xyl1-3' 'lad' 'lxr' 'xyl2' 'Seq1183' 'Seq4936' 'gal10' 'gal10-2' 'gal1' 'gal1-2' 'gal7'};

[a,b] = ismember(lacGenes,dataTable.genes);
lacMatrix = pearson_mat(b,b);
newVect = reshape(pearson_mat,nGenes*nGenes,1);
heatmap_obj = HeatMap(lacMatrix,'ColumnLabels',lacIDs,'RowLabels',lacIDs,'Colormap',redbluecmap,'Annotate',true);
h = plot((heatmap_obj));
h.FontSize = 18;
colormap(h,"bone")

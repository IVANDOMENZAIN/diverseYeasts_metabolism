leloir = [];%zeros(numel(csoure)*numel(condition),1);
oxredu = [];%zeros(numel(csoure)*numel(condition),1);
cd_str = {};
both   = [];
secrns = [];
sec_le = [];
sec_ox = [];
nadhUs = [];
nadpUs = [];
bothNAD = [];
sec_ndh = [];
sec_nph = [];
ethanol = [];
galactitol = [];
sorbose = [];
sorbitol = [];
tagatose = [];
j=1;
corrVectors = zeros(numel(model.mets),6);
pVals = zeros(numel(model.mets),6)+1;
for csoure = {'gal' 'lac'}
    for condition = {'clim' 'oxlim' 'nlim'}
        t = readtable(['../results/randomSampling_WT_' csoure{1} '_growth_' condition{1} '_solutions.txt'],'delimiter','\t');
        t = t(:,5:end);
        t = table2array(t);
        load('../models/candida_intermedia/cintGEM_curated.mat')
        [matrix] = getMetTurnOver(t,model);
        %matrix = t;
        %[coeff,score,latent,tsquared,explained] = pca(matrix');

        x = find(strcmp(model.metNames,'alpha-D-galactose 1-phosphate'));
        vector = matrix(x,:);%/max(matrix(x,:));
        idx1 = find(vector>0);
        leloir = [leloir; numel(idx1)];

        x      = find(strcmp(model.metNames,'galactitol'));
        x      = x(1);
        vector = matrix(x,:);%/max(matrix(x,:));
        idx2   = find(vector>0);
        oxredu = [oxredu; numel(idx2)];
        idx3   = intersect(idx1,idx2);
        both   = [both;numel(idx3)];
        cd_str = [cd_str; {[csoure{1} '_' condition{1}]}];

        u = find(strcmp(model.metNames,'galactitol'));
        u = u(model.metComps(u) == 3);
        vector = matrix(u,:);%/max(matrix(u,:));
        u = find(vector>0);
        galactitol = [galactitol;numel(u)];

        y = find(strcmp(model.metNames,'D-tagatose'));
        y = y(model.metComps(y) == 3);
        vector = matrix(y,:);%/max(matrix(y,:));
        y = find(vector>0);
        tagatose = [tagatose;numel(y)];

        z = find(strcmp(model.metNames,'L-sorbose'));
        z = z(model.metComps(z) == 3);
        vector = matrix(z,:);%/max(matrix(z,:));
        z = find(vector>0);
        sorbose = [sorbose;numel(z)];

        w = find(strcmp(model.metNames,'D-glucitol'));
        w = w(model.metComps(w) == 3);
        vector = matrix(w,:);%/max(matrix(w,:));
        w = find(vector>0);
        sorbitol = [sorbitol;numel(w)];

        v = find(strcmp(model.metNames,'ethanol'));
        v = v(model.metComps(v) == 3);
        vector = matrix(v,:);%/max(matrix(v,:));
        v = find(vector>0);
        ethanol = [ethanol;numel(v)];


        secretions = unique([u,y,z,w]);
        intLe = intersect(secretions,idx1);
        intOx = intersect(secretions,idx2);
        
        secrns = [secrns;numel(secretions)];
        sec_le = [sec_le;numel(intLe)];
        sec_ox = [sec_ox;numel(intOx)];

        x = find(strcmp(model.rxns,'ald_red_NADPH'));
        NADPH = find(t(x,:)~=0);
        y = find(strcmp(model.rxns,'ald_red_NADH'));
        NADH = find(t(y,:)~=0);
        z = intersect(NADPH,NADH);

        nadhUs = [nadhUs;numel(NADH)];
        nadpUs = [nadpUs;numel(NADPH)];
        bothNAD = [bothNAD;numel(z)];

        sec_nadh = intersect(secretions,NADH);
        sec_naph = intersect(secretions,NADPH);

        sec_ndh = [sec_ndh;numel(sec_nadh)];
        sec_nph = [sec_nph;numel(sec_naph)];
        
        NDPHMAT = t(:,NADPH);

        vector = zeros(numel(model.rxns),1)+1;
        for i=1:length(model.mets)
            %if ismember(i,NADPH)
                [R,P] = corrcoef(t(x,:),matrix(i,:));
                corrVectors(i,j) = R(1,2);
                pVals(i,j) = P(2,1);
            %end
        end
        topUp{j} = find(corrVectors(:,j)>0);
        topDn{j} = find(corrVectors(:,j)<-0);
        significant{j} = find(pVals(:,j)<=1E-1E-3);
        j = j+1;

    end
end
results = table(cd_str,leloir,oxredu,both,secrns,sec_le,sec_ox,nadhUs,nadpUs,bothNAD,sec_ndh,sec_nph);
secretion_summary = table(cd_str,ethanol,galactitol,tagatose,sorbitol,sorbose);

mets_gal_oxlim = sortrows(mets_gal_oxlim,'Var3','descend');
metsPos{3} = find(pVals(:,3)<=1E-3 & abs(corrVectors(:,3))>=0.1);
mets_gal_nlim= table(model.metNames(metsPos{3}),model.compNames(model.metComps(metsPos{3})),corrVectors(metsPos{3},3));
mets_gal_nlim = sortrows(mets_gal_nlim,'Var3','descend');
metsPos{4} = find(pVals(:,4)<=1E-3 & abs(corrVectors(:,4))>=0.1);
mets_lac_clim = table(model.metNames(metsPos{4}),model.compNames(model.metComps(metsPos{4})),corrVectors(metsPos{4},4));
mets_lac_clim = sortrows(mets_lac_clim,'Var3','descend');
metsPos{5} = find(pVals(:,5)<=1E-3 & abs(corrVectors(:,5))>=0.1);
mets_lac_oxlim = table(model.metNames(metsPos{5}),model.compNames(model.metComps(metsPos{5})),corrVectors(metsPos{5},5));
mets_lac_oxlim = sortrows(mets_lac_oxlim,'Var3','descend');
metsPos{6} = find(pVals(:,6)<=1E-3 & abs(corrVectors(:,6))>=0.1);
mets_lac_nlim = table(model.metNames(metsPos{6}),model.compNames(model.metComps(metsPos{6})),corrVectors(metsPos{6},6));
mets_lac_nlim = sortrows(mets_lac_nlim,'Var3','descend');
metsPos{3} = find(pVals(:,3)<=1E-3 & abs(corrVectors(:,3))>=0);
mets_gal_nlim= table(model.metNames(metsPos{3}),model.compNames(model.metComps(metsPos{3})),corrVectors(metsPos{3},3));
mets_gal_nlim = sortrows(mets_gal_nlim,'Var3','descend');
metsPos{5} = find(pVals(:,5)<=1E-3 & abs(corrVectors(:,5))>=0);
mets_lac_oxlim = table(model.metNames(metsPos{5}),model.compNames(model.metComps(metsPos{5})),corrVectors(metsPos{5},5));
mets_lac_oxlim = sortrows(mets_lac_oxlim,'Var3','descend');
metsPos{3} = find(pVals(:,3)<=1E-3 & abs(corrVectors(:,3))>=0);
mets_gal_nlim= table(model.metNames(metsPos{3}),model.compNames(model.metComps(metsPos{3})),corrVectors(metsPos{3},3));
metsPos{4} = find(pVals(:,4)<=1E-3 & abs(corrVectors(:,4))>=0.1);
mets_lac_clim = table(model.metNames(metsPos{4}),model.compNames(model.metComps(metsPos{4})),corrVectors(metsPos{4},4));
mets_lac_clim = sortrows(mets_lac_clim,'Var3','descend');
mets_gal_clim = sortrows(mets_gal_clim,'Var3','descend');
metsPos{5} = find(pVals(:,5)<=1E-3 & abs(corrVectors(:,5))>=0);
mets_lac_oxlim = table(model.metNames(metsPos{5}),model.compNames(model.metComps(metsPos{5})),corrVectors(metsPos{5},5));
mets_lac_oxlim = sortrows(mets_lac_oxlim,'Var3','descend');







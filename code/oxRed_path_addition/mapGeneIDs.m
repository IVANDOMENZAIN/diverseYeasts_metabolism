current = pwd;
%load model
load('../../models/candida_intermedia/cintGEM_oxido.mat');
%correct gene IDs (shorter strings)
model.genes = strrep(model.genes,'Candida_intermedia@','');
%We've run orthofinder using the FASTA files that were used for the model
%generation and the one for the RNAseq analysis (assembly steps)  to find
%orthologues between genes in the model (IDs: Seq_XXXX) and those in the
%RNAseq datasets (IDs: SGZXXXXXX.1)

%Let's check those results and map the genes in the model to the genes in
%the FASTA file. We'll store these new IDs in a new model field:
%model.orthologues

%checking the presence of the unassigned genes in the orthogroups
orthogroups     = readtable('../../orthoFinder/OrthoFinder/dataSEQs_vs_modelSEQs/Orthogroups/Orthogroups.txt','delimiter','\t');
[presence,idxs] = ismember(model.genes,orthogroups.model_Cint);
% it works!
idxs2 = find(presence);
idxs = idxs(presence);
model.orthologues = model.genes;
model.orthologues(idxs2) = orthogroups.data_Cint(idxs);
%overwrite the model
save('../../models/candida_intermedia/cintGEM_oxido.mat','model');

%In the FASTA file used for the RNAseq dataset there seem to be two kind of
%IDs for each gene (IDs: SGZXXXX and IDs CIC11Cxxxx), When checking the RNA
%dataset it seems that all genes are reported with the latter kind, so
%let's correct that in our model


%Open fasta file (the one used for RNAse1)
dataset = readtable('../../orthoFinder/data_Cint.txt','HeaderLines',0);
%Ignore lines with sequences
dataset = dataset(contains(dataset.ThisIsAFakeHeader,'>SGZ'),:);
%Get rid of the unnecessary characters in each column
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,'>','');
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,' [[Candida] intermedia]','');
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,' (mitochondrion)','');
%separate data by columns
column1 = [];
column2 = [];
for i=1:height(dataset)
    rowCell = strsplit(dataset.ThisIsAFakeHeader{i},' ');
    if length(rowCell)~=2
        warning('problem found in dataset.ThisIsAFakeHeader{i}')
        pause
    end
    column1 = [column1;rowCell(1)];
    column2 = [column2;rowCell(2)];
end
newDataset = table(column1,column2,'VariableNames',{'IDs_1' 'IDs_2'});
%Correct second column of IDs
newDataset.IDs_2 = strrep(newDataset.IDs_2,'CIC11C','CIC11T');
%Let's correct orthologues IDs in the model 
load('../../models/candida_intermedia/cintGEM_oxido.mat')
[presence,iB] = ismember(model.orthologues,newDataset.IDs_1);
iA = find(presence);
iB = iB(iB~=0);
%substitute IDs in the model
model.orthologues(iA) = newDataset.IDs_2(iB);
%There are some genes in the model that map to more than one sequence in
%the data, let's find them and substitute them
for i=1:length(model.orthologues)
    orthoID = model.orthologues{i};
    components = strsplit(orthoID,', ');
    if length(components)>1
        for j = 1:length(components)
            [~,iB] = ismember(components{j},newDataset.IDs_1);
            if iB>0
                components(j) = newDataset.IDs_2(iB);
            end
        end
    end
    model.orthologues{i} = strjoin(components,', ');
end

%When looking at the orthofinder results we've realized that there are some
%of the manually introduced genes (oxido-reductive pathway) that were
%already part of the model, let's correct this in the model and assign them
%with the correct ortholog ID

%LEt's substitute the IDs for the manually introduced genes
pos = find(strcmpi(model.orthologues,'xyl1'));
model.orthologues{pos} = 'CIC11T00000000334';
pos = find(strcmpi(model.orthologues,'xyl1_2'));
model.orthologues{pos} = 'CIC11T00000000893';
pos = find(strcmpi(model.orthologues,'xyl1_3'));
model.orthologues{pos} = 'CIC11T00000005922';

%WARNING: WE cannot find any orthologue fo lxr4 (trichoderma reesei)
%in the available sequence files for C. intermedia

%But let's save the model
save('../../models/candida_intermedia/cintGEM_oxido.mat','model')
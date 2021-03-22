%Open fasta file (the one used for RNAse1)
dataset = readtable('../orthoFinder/data_Cint.txt','HeaderLines',0);
%Ignore lines with sequences
dataset = dataset(contains(dataset.ThisIsAFakeHeader,'>SGZ'),:);
%Get rid of the unnecessary characters in each column
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,'>','');
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,' [[Candida] intermedia]','');
dataset.ThisIsAFakeHeader = strrep(dataset.ThisIsAFakeHeader,' (mitochondrion)','');

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
load('../models/candida_intermedia/cintGEM_oxido.mat')
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
save('../models/candida_intermedia/cintGEM_oxido.mat','model')
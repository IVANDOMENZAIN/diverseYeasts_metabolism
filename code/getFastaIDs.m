function newDataset = getFastaIDs
dataset = readtable('../orthoFinder/data_Cint.txt','HeaderLines',0);
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
end
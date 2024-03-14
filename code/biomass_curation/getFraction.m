function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid backbone');

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
if contains(rxnName,'lipid')
    subs = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F    = -sum(model.S(subs,fractionPos));   %g/gDW
else
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),model.mets{i});
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*(comps{pos,2}-18)/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

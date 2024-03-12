function [F,X] = getAApercentage(model,comps,X)
%Add up fraction:
fractionPos = find(strcmp(model.rxnNames,'protein pseudoreaction'));
comps = comps(strcmp(comps(:,3),'P'),:);
F = 0;
%Add up all components:
for i = 1:length(comps(:,1))
    pos = find(strcmp(model.mets,comps(i,1)));
    disp(model.metNames{pos})
    if numel(pos) == 1
        abundance = -model.S(pos,fractionPos)*(comps{i,2}-18)/1000;
        F         = F + abundance;
    end
    disp(num2str(abundance))
end
X = X + F;

perc = 0;
for i = 1:length(comps(:,1))
    pos = find(strcmp(model.mets,comps(i,1)));
    disp(model.metNames{pos})
    if numel(pos) == 1
        abundance = -model.S(pos,fractionPos)*(comps{i,2}-18)/1000;
    end
    percAA = abundance*100/X;
    disp(num2str(percAA))
    perc = percAA + perc;
end
end

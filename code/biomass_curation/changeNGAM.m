function modelMod =changeNGAM(modelMod,NGAM)
%Change GAM:
x = find(strcmpi(modelMod.rxnNames,'non-growth associated maintenance reaction'));
modelMod.lb(x) = NGAM;
end

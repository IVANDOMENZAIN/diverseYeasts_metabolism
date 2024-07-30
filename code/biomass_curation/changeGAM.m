function modelMod =changeGAM(modelMod,GAM)
%Change GAM:
xr_pos = strcmp(modelMod.rxns,'r_4041');
for i = 1:length(modelMod.mets)
    S_ix  = modelMod.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},modelMod.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        GAMpol = 0;
        modelMod.S(i,xr_pos) = sign(S_ix)*(GAM + GAMpol);
    end
end

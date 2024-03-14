function model = changeGAM(model,parameters,GAM)
%Change GAM:
xr_pos = strcmp(model.rxns,parameters.bioRxn);
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        GAMpol = 0;
        if isfield(parameters,'pol_cost')
            cost   = parameters.pol_cost;
            GAMpol = Ptot*cost(1) + Ctot*cost(2) + R*cost(3) + D*cost(4);
        end
        model.S(i,xr_pos) = sign(S_ix)*(GAM + GAMpol);
    end
end
end
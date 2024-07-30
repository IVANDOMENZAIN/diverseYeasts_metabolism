function model = rescalePseudoReaction(model,metName,f)
rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
if sum(rxnPos) == 1
    for i = 1:length(model.mets)
        S_ir   = model.S(i,rxnPos);
        isProd = strcmp(model.metNames{i},metName);
        if S_ir ~= 0 && ~isProd
            model.S(i,rxnPos) = f*S_ir;
        end
    end
end
end
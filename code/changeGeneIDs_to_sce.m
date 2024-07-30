function model = substituteGENEids(model)
for i=1:length(model.proteins)
    gene     = model.genes{i};
    ortholog = model.proteins{i};
    if ~isempty(ortholog)
        model.genes(i) = {ortholog};
        geneRxns = find(contains(model.grRules,gene));
        for j=1:length(geneRxns)
            model.grRules(geneRxns(j)) = strrep(model.grRules(geneRxns(j)),gene,ortholog);
        end
    end
    
end
end
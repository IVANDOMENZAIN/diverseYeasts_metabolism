load('../../models/candida_intermedia/cintGEM_oxido_curated.mat')
exp_val(1) = 0.66; %glucose
exp_val(2) = 0.59;
%unconstrain NGAM
x = find(strcmpi(model.rxnNames,'non-growth associated maintenance reaction'));
model.lb(x) = 0;
model.ub(x) = 1000;
%simulate growth on carbon-limited conditions (glucose)
model = changeMedia_batch(model,'D-glucose exchange',1);
bio_pos = find(strcmp(model.rxnNames,'growth'));
glu_pos = find(strcmp(model.rxnNames,'D-glucose exchange'));
model = setParam(model,'obj',bio_pos,1);
sol_glu   = solveLP(model,1);
bioYield_glu = sol_glu.x(bio_pos)/abs(sol_glu.x(glu_pos)*0.18);
%simulate growth on carbon-limited conditions (lactose)
model = changeMedia_batch(model,'lactose exchange',1);
lac_pos = find(strcmp(model.rxnNames,'lactose exchange'));
model = setParam(model,'obj',bio_pos,1);
sol   = solveLP(model,1);
bioYield_lac = sol.x(bio_pos)/abs(sol.x(lac_pos)*0.3429);
error(1) = (bioYield_glu-exp_val(1))/exp_val(1);
error(2) = (bioYield_glu-exp_val(2))/exp_val(2);
%it seems that the biomass yield on glucose is underpredicted (5%) and
%lactose overpredicted by 6%. let's check the stoichiometries in Leloir
%pathway
leloirGenes = {'YBR019C' 'YBR020W' 'YBR018C' 'YHL012W' 'YKL127W' 'YCL040W'};
leloirShort = {'GAL10' 'GAL1' 'GAL7' 'YHL012W' 'PGM1' 'GLK1'};

%Check the presence of leloir's genes in cintGEM, print rxn formulas
notPresent = [];
for i=1:length(leloirGenes)
    idx = find(strcmpi(model.proteins,leloirGenes{i}));
    if ~isempty(idx)
        for j=1:length(idx)
            disp(leloirShort{i})
            disp(model.genes(idx(j)))
            rxnIdxs = find(model.rxnGeneMat(:,idx(j)));
            printModel(model,rxnIdxs)
        end
    else
        disp(['Gene: ' leloirGenes{i} ' is not present in cIntGEM']) 
        notPresent = [notPresent;leloirGenes(i)];
    end
end
%we found a weird reaction (probably from cerevisiae) for GAL7, 'r_0459',
%let's check if this offers an energy advantage that might have an impact
%on bio yield on lactose
modelConst = setParam(model,'lb','r_0459',0);
modelConst = setParam(modelConst,'ub','r_0459',0);
modelConst = setParam(modelConst,'lb','r_4527',0);
modelConst = setParam(modelConst,'ub','r_4527',0);

modelConst = changeMedia_batch(modelConst,'lactose exchange',1);
lac_pos = find(strcmp(modelConst.rxnNames,'lactose exchange'));
modelConst = setParam(modelConst,'obj',bio_pos,1);
sol_lac   = solveLP(modelConst,1);
bioYield_lac = sol_lac.x(bio_pos)/abs(sol_lac.x(lac_pos)*0.3429);
error(1) = (bioYield_glu-exp_val(1))/exp_val(1);
error(2) = (bioYield_glu-exp_val(2))/exp_val(2);
%compare flux dist
formulas = constructEquations(model);
FCs = (sol_lac.x+1E-6)./(sol_glu.x+1E-6);
fluxTable = table(model.rxns,model.rxnNames,formulas,sol_glu.x,sol_lac.x,FCs,model.grRules);
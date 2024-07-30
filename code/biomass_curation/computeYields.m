load('../../models/candida_intermedia/cintGEM_curated.mat')
exp_val(1) = 0.66; %glucose
exp_val(2) = 0.59;
%unconstrain NGAM
%x = find(strcmpi(model.rxnNames,'non-growth associated maintenance reaction'));
%model.lb(x) = 0;
%model.ub(x) = 1000;
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
%block cytosolic ATPase
%modelConst = setParam(modelConst,'ub','r_0227',0);
%modelConst = setParam(modelConst,'lb','r_0227',0);

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
fluxTable_glcVsLac = table(model.rxns,model.rxnNames,formulas,sol_glu.x,sol_lac.x,FCs,model.grRules);
galacPos = find(strcmpi(model.metNames,'galactitol'));
disp(model.metNames(galacPos))
disp(model.metComps(galacPos))
galacRxn = find(model.S(galacPos(2),:));
constructEquations(model,galacRxn)
%simulate max. galactitol production
modelConst = setParam(modelConst,'obj',galacPos,1);
modelConst.ub(galacRxn) =1000;
sol_lac   = solveLP(modelConst,1);
galactitol = sol_lac.x(galacRxn)/abs(sol_lac.x(lac_pos));

%simulate galactitol accumulation
modelConst = setParam(modelConst,'obj',galacRxn,1);
sol_lac   = solveLP(modelConst,1);
galactitol = sol_lac.x(galacRxn)/abs(sol_lac.x(lac_pos));
printFluxes(modelConst,sol_lac.x,true)
%simulate growth of GAL mutant
GALgenes = {'Seq_1935' 'Seq_4294' ... %gal1
            'Seq_3460' ... %gal10
            'Seq_2479' 'Seq_3332'};
GALmutant = removeGenes(modelConst,GALgenes,true,false,true);
%simulate growth
GALmutant = setParam(GALmutant,'obj',bio_pos,1);
sol_lac_mut   = solveLP(GALmutant,1);

modelConst = setParam(modelConst,'obj',bio_pos,1);
sol_lac_WT   = solveLP(modelConst,1);
indxs = (abs(sol_lac_mut.x)+abs(sol_lac_WT.x))>0;
formulas = constructEquations(modelConst);
FCs = (sol_lac_mut.x+1E-6)/(sol_lac_WT.x+1E-6);
fluxTable = table(modelConst.rxns,modelConst.rxnNames,modelConst.grRules,formulas,sol_lac_WT.x,sol_lac_mut.x);
%simulate growth on galactitol
GALmutant = setParam(GALmutant,'obj',bio_pos,1);
GALmutant = changeMedia_batch(GALmutant,'galactitol exchange',1);
sol_lac_mut   = solveLP(GALmutant,1);
modelConst = setParam(modelConst,'obj',bio_pos,1);
modelConst = changeMedia_batch(modelConst,'galactitol exchange',1);
sol_lac_WT   = solveLP(modelConst,1);
indxs = (abs(sol_lac_mut.x)+abs(sol_lac_WT.x))>0;
FCs = (sol_lac_mut.x+1E-6)/(sol_lac_WT.x+1E-6);
fluxTable = table(modelConst.rxns(indxs),modelConst.rxnNames(indxs),modelConst.grRules(indxs),formulas(indxs),sol_lac_WT.x(indxs),sol_lac_mut.x(indxs),FCs(indxs));

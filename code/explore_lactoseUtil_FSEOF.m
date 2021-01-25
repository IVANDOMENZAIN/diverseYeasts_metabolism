%explore_CintGEM_wFSEOF
load('../models/candida_intermedia/cintGEM_oxido.mat')
%identify relevant rxns
biomassRxn = find(strcmpi(model.rxnNames,'biomass pseudoreaction'));
biomassRxn = model.rxns(biomassRxn);
%Rxn ID for L-glycine exchange reaction (high production capabilities)
targetRxn  = 'r_1810';
%retrieve lactose exchange reaction
lacRxn = find(strcmpi(model.rxnNames,'lactose exchange'));
lacRxn = model.rxns{lacRxn};
%Set media constraints with lactose as a carbon source
modelLac = changeMedia(model,lacRxn,1);
%Verify growth 
modelLac = setParam(modelLac,'obj',biomassRxn,1);
solution = solveLP(modelLac,1);
printFluxes(modelLac,solution.x)
%Everything's fine so far

%Let's try to use the FSEOF function
mkdir('../results/FSEOF_cint/')
outputFile = '../results/FSEOF_cint/targets_test.txt';
%targets=FSEOF(modelLac,biomassRxn,targetRxn,10,0.9,outputFile);
%It works! (for production L-glycine)
%LEt's use our function to study fdifferential usage of lactose utilization
%pathways subject to a fixed biomass production

%Fix biomass
modelLac = setParam(modelLac,'lb',biomassRxn,0.9999*0.1861);
modelLac = setParam(modelLac,'ub',biomassRxn,1.0001*0.1861);
%Block alternative reaction to galactokinase
modelLac = setParam(modelLac,'lb','r_4222',0);
modelLac = setParam(modelLac,'ub','r_4222',0);

galactokinase = 'r_0458';
oxRed_path = 'xyl_hex_red';
outputFile = '../results/FSEOF_cint/upReg_fluxes_glcK_2_OxRed.txt';
%targets=FSEOF(modelLac,galactokinase,oxRed_path,10,0.99,outputFile);

%Block ald_red_NADH to force usage of the NADPH-specific isoform
modelLac = setParam(modelLac,'lb','ald_red_NADH',0);
modelLac = setParam(modelLac,'ub','ald_red_NADH',0);

outputFile = '../results/FSEOF_cint/upReg_fluxes_glcK_2_OxRed_NADPH.txt';
targets    = extendedFSEOF(modelLac,galactokinase,oxRed_path,10,0.99,outputFile);


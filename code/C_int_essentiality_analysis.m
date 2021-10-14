%Essentiality analysis
load('..\models\candida_intermedia\cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x) %It is growing on lactose!
%Making sure...
lacRxn = find(strcmpi(model.rxnNames,'lactose exchange'));
lacRxn = model.rxns{lacRxn};


%GAL1_2 deletion doesn't grow on lactose IRL
%Let's block it
galactokinase = find(strcmp(model.rxns,'r_4222'))
model = setParam(model,'lb','r_4222',0);
model = setParam(model,'ub','r_4222',0);
%apply the changes
sol = solveLP(model,1);
printModel(model,galactokinase) %galactokinase is blocked
printFluxes(model,sol.x) %but it still grows
%THAT IS A FALSE POSITIVE

%Novel cluster - GAL1_2, GAL10_2, Lac9 and XYL1_2, 
%all together when deleted, resulted in no growth on lactose and also galactose.
model = setParam(model,'lb','r_4222',0);
model = setParam(model,'ub','r_4222',0);


%We still need to introduce the doubled isoenzymes that are in the
%oxidorreductive pathway, cause they are still not there.

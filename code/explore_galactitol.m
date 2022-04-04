%explore galactitol production in C. intermedia model
load('..\models\candida_intermedia\cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x)
%Model is already growing on lactose, growth is 0.20682 if lb of lac_ex = 1

%it looks like other yeasts (such as galactose-fermenting S.
%cerevisiae strains) display intracellular accumulation of galactitol (https://doi.org/10.1002/bit.21890), so

%Let's see if galactitol is produced
gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
sol.x(gtoh_x)
%No galactitol produced at all

%Let's just see if we have flux along the redox pathway
rxns = {'ald_red_NADPH' 'r_4983' 'xyl_hex_red' 'r_5174' 'r_0323'};
fluxes=haveFlux(model,1E-6,rxns);
% They do

%let's see what happens when we block the leloir pathway:
index = find(contains(model.rxns,'r_0458')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 0;
index = find(contains(model.rxns,'r_4222')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 0;
sol = solveLP(model,1);
printFluxes(model,sol.x) %Growth: 0.19524

%Let's see if galactitol is produced now
gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
sol.x(gtoh_x)
%It was! Flux was 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now, let's just reset everything and increase the consumption of lactose
load('..\models\candida_intermedia\cintGEM_oxido_AVR.mat')
%Let's see if everything was reset
sol = solveLP(model,1);
sol.x(index) %Leloir unblocked
sol.x(gtoh_x) %galactitol back to no production.

% If galactitol is an overflow metabolite, high fluxes of lactose should
% induce production.
sol.x(index) %Leloir is still unblocked
glu_ex = find(contains(model.rxnNames,'glucose exchange'));
printModel(model,glu_ex) %We don't have glucose exchange, only lactose:
lacEx = find(strcmpi(model.rxns,'lac_ex'));
printModel(model,lacEx) % lb is at -1, let's see what happens if we increase it to -1000
model.lb(lacEx)   = -1000;
sol = solveLP(model,1);
sol.x(gtoh_x) %We still have no production at high lactose rates

%Let's block leloir again
%let's see what happens when we block the leloir pathway:
index = find(contains(model.rxns,'r_0458')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 0;
index = find(contains(model.rxns,'r_4222')); %Galactokinase
model.lb(index)  = 0;
model.ub(index)  = 0;
sol = solveLP(model,1);
printFluxes(model,sol.x) %Growth: 0.19524
haveFlux(model,1E-6,'r_0458')
haveFlux(model,1E-6,'r_4222') %no flux through leloir anymore

%Let's look at the import of lactose now
printModel(model,lacEx) %bounds are [-1 1000]
sol.x(gtoh_x) %Galactitol production is 1
% let's now open lactose bounds a little bit
model.lb(lacEx)   = -500;
model.lb(glu_ex)   = 0;
model.ub(glu_ex)   = 0;
sol = solveLP(model,1);
printFluxes(model,sol.x) %Lactose absorption is 500
sol.x(gtoh_x) %Galactitol too
%Let's open up even more
model.lb(lacEx)   = -1000;
sol = solveLP(model,1);
printFluxes(model,sol.x) %Lactose absorption is 670.91
sol.x(gtoh_x) %Galactitol too

%Interesting. Let's see if we can model lactose intake, growth and galactitol
%production as the cell takes more and more lactose

% Just looking for growth reaction:
grow = find(contains(model.rxns,'r_2111'));

%Let's check glycolysis
hexo = find(contains(model.rxnNames,'D-fructose:ATP'));
sol.x(hexo)

sorbitol = find(contains(model.rxnNames,'sorbitol'));
sol.x(sorbitol)

%find out where we're getting out sorbitol from



%what happens to galactose as soon as it enters the cell

model.rxnNames(find(contains(model.rxnNames,'galactosidase')))














hexulose = find(contains(model.rxnNames,'arabitol'));

%Let's just shut down the hexokinase for fructose, which we think is being
%used by the OR pathway
model.lb(lacEx)   = 0;
model.ub(lacEx)   = 0;
%WHERE THE HELL IS GLUCOSE GOING

f6p = find(contains(model.rxnNames,'phosphofructokinase'));
sol.x(f6p)


pep = find(contains(model.rxnNames,'enolase'));
sol.x(pep)

pyruvate = find(contains(model.rxnNames,'pyruvate kinase'));
sol.x(pyruvate)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's go back to echange of 1 on lactose
model.lb(lacEx) = -1;
model.lb(glu_ex) = 0;
sol = solveLP(model,1);
printFluxes(model,sol.x)

%We see that after a certain point growth is slowed down and so is the
%intake of lactose. (even though we're not actually doing anything on C.
%intermedia, we can try looking into this now and then doing that stuff on
%the improved model.










% Kamesh mentioned that peaks may overlap and that maybe galactitol wasn't
% what appeared on the spectrum, but sorbitol or another sugar alcohol

%Let's find out if any of the products in the pathway is secreted
model.rxnNames(find(contains(model.rxnNames,'D-fructose exchange')))
%We have a fructose exchange, but it's not a sugar alcohol.
%What happens to everything if we add a galactitol exchange reaction?

% Define reaction equation
exchange_gtOH = 'galactitol[c] => '
rxnsToAdd.equations = {exchange_gtOH}; 
% Define reaction names
rxnsToAdd.rxns     = {'gtOH_ex'};
rxnsToAdd.rxnNames = {'galactitol exchange'};
% Define objective and bounds
rxnsToAdd.c  = [0];
rxnsToAdd.lb = [0];
rxnsToAdd.ub = [1000];

%Add to model
model_gtOH = addRxns(model,rxnsToAdd,3);

%Evaluate if rxn can carry flux
I = haveFlux(model_gtOH,1E-6,'gtOH_ex') %It can!!! Let's take it for a 
%test run with high lactose

model_gtOH.lb(lacEx) = -1;
sol = solveLP(model_gtOH,1);
gtol = find(contains(model_gtOH.rxns,'gtOH_ex'));
grow = find(contains(model_gtOH.rxns,'r_2111'));
model_gtOH.lb(gtol)  = 0;
model_gtOH.ub(gtol)  = 1000;

index = find(contains(model_gtOH.rxns,'r_0458')); %Galactokinase
model_gtOH.lb(index)  = 0;
model_gtOH.ub(index)  = 0;
index = find(contains(model_gtOH.rxns,'r_4222')); %Galactokinase
model_gtOH.lb(index)  = 0;
model_gtOH.ub(index)  = 0;
sol = solveLP(model_gtOH,1);
printFluxes(model_gtOH,sol.x) %Growth: 0.19524

sol = solveLP(model_gtOH,1);
printFluxes(model_gtOH,sol.x) %It grows the same as without exch rxn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_growth=setParam(model_gtOH,'obj',{'r_2111'},1);
sol=solveLP(model_growth,1);
printFluxes(model_growth, sol.x, true, 10^-7);

model_gtol=setParam(model_gtOH,'obj',{'gtOH_ex'},1);
solgtol=solveLP(model_gtol,1);
printFluxes(model_gtol, solgtol.x, true, 10^-7);


followChanged(model_gtol,sol.x,solgtol.x, 50, 0.5, 0.5);

load 'pcPathway.mat' pathway;

drawMap('Growth vs galactitol',pathway,model_gtol,solgtol.x,sol.x,model_growth,'test2.pdf',10^-5);



% 
% printModelStats(model_gtOH,true,true)
% 
% %D-glucose appears as a dead-end metabolite
% %The reactions we introduced are not elementally balanced
% 
% % model=simplifyModel(model,true,false,true,true);
% % printModelStats(model,true,true)
% 
% printFluxes(model_gtOH, sol.x, false, 10^-7);
% 
% sol=solveLP(model,1);
% printFluxes(model, sol.x, false, 10^-7);
% 
% model=setParam(model,'obj',{'r_2111'},1);
% sol=solveLP(model,1);
% printFluxes(model, sol.x, true, 10^-7);





%IF FLUXES ARE CMOLES THEN WHY ARE WE USING 1 LACTOSE AND NOT 2

lactose_grad = [];
growth_grad = [];
gtOH_grad = [];

for i = 1:1000
    model_gtOH.lb(lacEx) = -i;
    sol = solveLP(model_gtOH,1);
    lactose_g = -sol.x(lacEx);
    growth_g = sol.x(grow);
    gtOH_g = sol.x(gtol);
    %Let's fill the matrix with our outputs
    lactose_grad = [lactose_grad; lactose_g];
    growth_grad = [growth_grad; growth_g];
    gtOH_grad = [gtOH_grad;gtOH_g];
    i
end

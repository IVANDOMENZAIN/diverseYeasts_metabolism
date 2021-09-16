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

% Galactitol is an overflow metabolite, so high fluxes of lactose should
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

%Define empty matrix to fill during loop, this will be used for plot
%later

lactose = [];
growth_lac = [];
gtOH_lac = [];
i2 = [];

for i = 1:1000
    model.lb(lacEx) = -i;
    sol = solveLP(model,1);
    lactose2 = -sol.x(lacEx);
    growth2 = sol.x(grow);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    lactose = [lactose; lactose2];
    growth_lac = [growth_lac; growth2];
    gtOH_lac = [gtOH_lac;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(1,3,1)
plot(i2, lactose)
xlabel("lactose bounds")
ylabel("lactose intake")
legend("lactose intake as we open the bounds")
subplot(1,3,3)
plot(i2, gtOH_lac)
xlabel("lactose bounds")
ylabel("gtOH production")
legend("gtOH production as we open the bounds")
subplot(1,3,2)
plot(i2, growth_lac)
xlabel("lactose bounds")
ylabel("growth")
legend("growth as we open the bounds")

%Let's now do glucose instead of lactose
sol.x(index) %Leloir is still blocked
glu_ex = find(contains(model.rxnNames,'glucose exchange'));
model.lb(glu_ex)   = -1000;
printModel(model,glu_ex) %We don't have glucose exchange, only lactose:
lacEx = find(strcmpi(model.rxns,'lac_ex'));
model.lb(lacEx)   = 0;
printModel(model,lacEx) % lb is at -1, let's see what happens if we increase it to -1000
%Let's solve
sol = solveLP(model,1);
sol.x(gtoh_x) %No galactitol again

%Define empty matrix to fill during loop, this will be used for plot
%later

glucose = [];
growth_glu = [];
gtOH_glu = [];
i2 = [];

for i = 1:1000
    model.lb(glu_ex) = -i;
    sol = solveLP(model,1);
    glucose2 = -sol.x(glu_ex);
    growth2 = sol.x(grow);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    glucose= [glucose; glucose2];
    growth_glu = [growth_glu; growth2];
    gtOH_glu = [gtOH_glu;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(1,3,1)
plot(i2, glucose)
xlabel("glucose bounds")
ylabel("glucose intake")
legend("glucose intake as we open the bounds")
subplot(1,3,3)
plot(i2, gtOH_glu)
xlabel("glucose bounds")
ylabel("gtOH production")
legend("gtOH production as we open the bounds")
subplot(1,3,2)
plot(i2, growth_glu)
xlabel("glucose bounds")
ylabel("growth")
legend("growth as we open the bounds")

%Now let's do it for the both of them!

glucose_2 = [];
lactose_2 = [];
growth_2 = [];
gtOH_2 = [];
i2 = [];

for i = 1:1000
    model.lb(glu_ex) = -i;
    model.lb(lacEx) = -i;
    sol = solveLP(model,1);
    glucose2 = -sol.x(glu_ex);
    lactose2 = -sol.x(lacEx);
    growth2 = sol.x(grow);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    glucose_2 = [glucose_2; glucose2];
    lactose_2 = [lactose_2; lactose2];
    growth_2 = [growth_2; growth2];
    gtOH_2 = [gtOH_2;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(2,2,1)
plot(i2, glucose_2)
xlabel("glucose bounds")
ylabel("glucose intake")
legend("glucose intake as we open both bounds")
subplot(2,2,2)
plot(i2, lactose_2)
xlabel("lactose bounds")
ylabel("lactose intake")
legend("lactose intake as we open both bounds")
subplot(2,2,3)
plot(i2, growth_2)
xlabel("bounds")
ylabel("growth")
legend("growth as we open the bounds")
subplot(2,2,4)
plot(i2, gtOH_2)
xlabel("bounds")
ylabel("gtOH production")
legend("gtOH production as we open the bounds")

%Now let's see some gradients

glucose_grad = [];
lactose_grad = [];
growth_grad = [];
gtOH_grad = [];

for i = 1:1000
    model.lb(glu_ex) = -i;
    model.lb(lacEx) = -(1000-i);
    sol = solveLP(model,1);
    glucose_g = -sol.x(glu_ex);
    lactose_g = -sol.x(lacEx);
    growth_g = sol.x(grow);
    gtOH_g = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    glucose_grad = [glucose_grad; glucose_g];
    lactose_grad = [lactose_grad; lactose_g];
    growth_grad = [growth_grad; growth_g];
    gtOH_grad = [gtOH_grad;gtOH_g];
    i
end

%Now let's generate plots to see the results
subplot(2,2,1)
plot(i2, glucose_grad)
xlabel("glucose bounds")
ylabel("glucose intake")
legend("glucose intake as we open glucose bounds and close lactose bounds")
subplot(2,2,2)
plot(i2, lactose_grad)
xlabel("lactose bounds")
ylabel("lactose intake")
legend("lactose intake as we open glucose bounds and close lactose bounds")
subplot(2,2,3)
plot(i2, growth_grad)
xlabel("bounds")
ylabel("growth")
legend("growth as we open glucose bounds and close lactose bounds")
subplot(2,2,4)
plot(i2, gtOH_grad)
xlabel("bounds")
ylabel("gtOH production")
legend("gtOH production as we open glucose bounds and close lactose bounds")


%Reverse gradient now

%Now let's see some gradients

glucose_darg = [];
lactose_darg = [];
growth_darg = [];
gtOH_darg = [];

for i = 1:1000
    model.lb(glu_ex) = -(1000-i);
    model.lb(lacEx) = -i;
    sol = solveLP(model,1);
    glucose_t = -sol.x(glu_ex);
    lactose_t = -sol.x(lacEx);
    growth_t = sol.x(grow);
    gtOH_t = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    glucose_darg = [glucose_darg; glucose_t];
    lactose_darg = [lactose_darg; lactose_t];
    growth_darg = [growth_darg; growth_t];
    gtOH_darg = [gtOH_darg;gtOH_t];
    i
end

%Now let's generate plots to see the results
subplot(2,2,1)
plot(i2, glucose_darg)
xlabel("glucose bounds")
ylabel("glucose intake")
legend("glucose intake as we open lactose bounds and close glucose bounds")
subplot(2,2,2)
plot(i2, lactose_darg)
xlabel("lactose bounds")
ylabel("lactose intake")
legend("lactose intake as we open lactose bounds and close glucose bounds")
subplot(2,2,3)
plot(i2, growth_darg)
xlabel("bounds")
ylabel("growth")
legend("growth as we open lactose bounds and close glucose bounds")
subplot(2,2,4)
plot(i2, gtOH_darg)
xlabel("bounds")
ylabel("gtOH production")
legend("gtOH production as we open lactose bounds and close glucose bounds")

%Let's now just plot growth comparisons
subplot(3,2,1)
plot(i2, growth_glu)
xlabel("glucose bounds")
ylabel("growth")
legend("growth as we open the bounds")
subplot(3,2,2)
plot(i2, growth_lac)
xlabel("lactose bounds")
ylabel("growth")
legend("growth as we open the bounds")
subplot(3,2,3)
plot(i2, growth_2)
xlabel("bounds")
ylabel("growth")
legend("growth as we open both bounds")
subplot(3,2,4)
plot(i2, growth_grad)
xlabel("bounds")
ylabel("growth")
legend("growth as we open glucose bounds and close lactose bounds")
subplot(3,2,5)
plot(i2, growth_darg)
xlabel("bounds")
ylabel("growth")
legend("growth as we open lactose bounds and close glucose bounds")



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

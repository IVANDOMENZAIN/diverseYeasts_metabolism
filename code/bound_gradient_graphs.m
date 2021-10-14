%Let's see some simulations
load('..\models\candida_intermedia\cintGEM_oxido_AVR.mat')
%See how exchange fluxes look
sol = solveLP(model,1);
printFluxes(model,sol.x)
%Model is already growing on lactose, growth is 0.20682 if lb of lac_ex = 1

lacRxn = find(strcmpi(model.rxnNames,'lactose exchange'));
lacRxn_i = model.rxns{lacRxn};
glu_ex = find(contains(model.rxnNames,'glucose exchange'));
glu_ex_i = model.rxns{glu_ex};
gtoh_x = find(contains(model.rxnNames,'aldose reductase (NADH)'));
gtoh_x_i = model.rxns{gtoh_x};
growth = find(contains(model.rxns,'r_2111'));
gal_ex = find(contains(model.rxnNames,'D-galactose exchange'));
gal_ex_i = model.rxns{gal_ex};


%Let's make sure the model is growing on lactose

model = changeMedia(model,lacRxn_i,15);
printModel(model,lacRxn)
printModel(model,glu_ex)
printModel(model,gtoh_x)
sol = solveLP(model,1);
printFluxes(model,sol.x)

model = changeMedia(model,glu_ex_i,15);
printModel(model,lacRxn)
printModel(model,glu_ex)
printModel(model,gtoh_x)
printFluxes(model,sol.x)



%Define empty matrix to fill during loop, this will be used for plot
%later

lactose = [];
growth_lac = [];
gtOH_lac = [];
i2 = [];

for i = 1:1000
    model.lb(lacRxn) = -i;
    sol = solveLP(model,1);
    lactose2 = -sol.x(lacRxn);
    growth2 = sol.x(growth);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    lactose = [lactose, lactose2];
    growth_lac = [growth_lac; growth2];
    gtOH_lac = [gtOH_lac;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(1,2,1)
plot(i2, growth_lac)
xlabel("lactose bounds")
ylabel("growth")
legend("growth as we open the bounds")
subplot(1,2,2)
plot(i2, gtOH_lac)
xlabel("lactose bounds")
ylabel("gtOH production")
legend("gtOH production as we open the bounds")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's now do glucose instead of lactose
sol.x(index) %Leloir is still blocked
glu_ex = find(contains(model.rxnNames,'glucose exchange'));
model.lb(glu_ex)   = -1000;
printModel(model,glu_ex) %We don't have glucose exchange, only lactose:
lacEx = find(strcmpi(model.rxns,'lac_ex'));
model.lb(lacRxn)   = 0;
printModel(model,lacEx) % lb is at -1, let's see what happens if we increase it to -1000
%Let's solve
sol = solveLP(model,1);
sol.x(gtoh_x) %No galactitol again

%Define empty matrix to fill during loop, this will be used for plot
%later

glucose = zeros(1000,1);
growth_glu = zeros(1000,1);
gtOH_glu = zeros(1000,1);
i2 = zeros(1000,1);

for i = 1:1000
    model.lb(glu_ex) = -i;
    sol = solveLP(model,1);
    glucose2 = -sol.x(glu_ex);
    growth2 = sol.x(growth);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    glucose= [glucose; glucose2/180.56];
    growth_glu = [growth_glu; growth2];
    gtOH_glu = [gtOH_glu;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(1,2,2)
plot(i2, gtOH_glu)
xlabel("glucose bounds")
ylabel("gtOH production rate")
legend("gtOH production as we open the bounds")
subplot(1,2,1)
plot(i2, growth_glu)
xlabel("glucose bounds")
ylabel("growth rate")
legend("growth as we open the bounds")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now with galactose

galactose = zeros(1000,1);
growth_gal = zeros(1000,1);
gtOH_gal = zeros(1000,1);
i2 = zeros(1000,1);

for i = 1:1000
    model.lb(gal_ex) = -i;
    sol = solveLP(model,1);
    galactose2 = -sol.x(gal_ex);
    growth2 = sol.x(growth);
    gtOH2 = sol.x(gtoh_x);
    %Let's fill the matrix with our outputs
    galactose = [galactose; galactose2];
    growth_gal = [growth_glu; growth2];
    gtOH_gal = [gtOH_glu;gtOH2];
    i2 = [i2;i];
    i
end

%Now let's generate plots to see the results
subplot(1,2,2)
plot(i2, gtOH_gal)
xlabel("galactose bounds")
ylabel("gtOH production rate")
legend("gtOH production as we open the bounds")
subplot(1,2,1)
plot(i2, growth_gal)
xlabel("galactose bounds")
ylabel("growth rate")
legend("growth as we open the bounds")





%Now let's do it for the both of them!

glucose_2 = zeros(1000,1);
lactose_2 = zeros(1000,1);
growth_2 = zeros(1000,1);
gtOH_2 = zeros(1000,1);
i2 = zeros(1000,1);

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
i2 = [];

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
    i2 = [i2;i];    
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
i2 = [];   
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
    i2 = [i2;i];   
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


load('../models/candida_intermedia/cintGEM_curated.mat')
nSamples = 10000;
flux_threshold = 1E-3;
%open growth set it as objective function, and set a unit lactose uptake
%rate
model = setParam(model,'lb','r_4041',0);
model = setParam(model,'ub','r_4041',1000);
model = setParam(model,'obj','r_4041',1);
model = changeMedia_batch(model,'lactose exchange',1);
x = find(strcmp(model.rxnNames,'lactose exchange'));
model.lb(x) = -1;
model.ub(x) = 0;
%get relevant reaction positions in the WT model
indexes(1) =  x;
indexes(2) =  find(strcmp(model.rxnNames,'oxygen exchange'));
indexes(3) =  find(strcmp(model.rxnNames,'growth'));
indexes(4) =  find(strcmp(model.rxnNames,'galactitol exchange'));
indexes(5) =  find(strcmp(model.rxnNames,'carbon dioxide exchange'));
indexes(6) =  find(strcmp(model.rxnNames,'ethanol exchange'));
indexes(7) =  find(strcmp(model.rxnNames,'D-tagatose exchange'));
indexes(8) =  find(strcmp(model.rxnNames,'L-sorbose exchange'));

%get GAL mutant model
mutant = getGALmutant(model);
solWT = solveLP(model,1);
solMT = solveLP(mutant,1);
%verify growth
maxG = solWT.x(find(model.c));
maxGM = solMT.x(find(mutant.c));
%get relevant reaction positions in the mutant model
indexesM(1) =  x;

indexesM(2) =  find(strcmp(mutant.rxnNames,'oxygen exchange'));
indexesM(3) =  find(strcmp(mutant.rxnNames,'growth'));
indexesM(4) =  find(strcmp(mutant.rxnNames,'galactitol exchange'));
indexesM(5) =  find(strcmp(mutant.rxnNames,'carbon dioxide exchange'));
indexesM(6) =  find(strcmp(mutant.rxnNames,'ethanol exchange'));
indexesM(7) =  find(strcmp(mutant.rxnNames,'D-tagatose exchange'));
indexesM(8) =  find(strcmp(mutant.rxnNames,'L-sorbose exchange'));

%initialize variables 
maxO2 = 1;
iterations = 100;
delta = maxO2/iterations;
model.lb(indexes(3)) = 0;
model.ub(indexes(3)) = 1000;
mutant.lb(indexesM(3)) = 0;
mutant.ub(indexesM(3)) = 1000;
exchangesMT = zeros(iterations,numel(indexes));
exchangesWT = zeros(iterations,numel(indexesM));
%iterate through increasing oxygen uptake levels and maximize growth on
%both models, store exchange fluxes
model = setParam(model,'ub','r_1883',0); 
mutant = setParam(mutant,'ub','r_1883',0); 

for i=1:iterations
    model = setParam(model,'obj','r_4041',1);
    mutant = setParam(mutant,'obj','r_4041',1);
    mutant = setParam(mutant,'obj','r_4041',1);
    %block ethanol production
    %model.ub(indexes(6)) = 0;
    %mutant.ub(indexesM(6)) = 0;

        %block butanediol production
    model = setParam(model,'ub','r_1549',0);
    mutant = setParam(mutant,'ub','r_1549',0);
    model = setParam(model,'ub','r_1875',0);
    mutant = setParam(mutant,'ub','r_1875',0);
    %set oxygen limitations
    model.lb(indexes(2)) = -(i)*maxO2/iterations;
    model.ub(indexes(2)) = 0;
    mutant.lb(indexesM(2)) = -(i)*maxO2/iterations;
    mutant.ub(indexesM(2)) = 0;
    %verify C source uptake
    model.lb(x) = -1;
    mutant.lb(x) = -1;
    model.ub(x) = -0.9999;
    mutant.ub(x) = -0.9999;
    %solve
    solWT = solveLP(model,1);
    solMT= solveLP(mutant,1);
    if ~isempty(solWT.x) & ~isempty(solMT.x)
        printFluxes(model,solWT.x,true)
        disp(' ')
        printFluxes(mutant,solMT.x,true)
        exchangesWT(i,:) = solWT.x(indexes)';
        exchangesMT(i,:) = solMT.x(indexesM)';
    end
end
%Visualize results for mutant model
matrix = abs(exchangesMT);
figure
plot(matrix(:,2),matrix(:,3),matrix(:,2),matrix(:,4),matrix(:,2),matrix(:,5),matrix(:,2),matrix(:,6),matrix(:,2),matrix(:,7),matrix(:,2),matrix(:,8),'LineWidth',4);
xlabel('Oxygen uptake [mmol/gDw h]','FontSize',18)
ylabel('Exchange flux [mmol/gDw h] - [1/h]','FontSize',18)
ylim([0 4])
legend(model.rxnNames(indexesM(3:end)))
%title('mutant galactose')
%saveas(fig,['../results/figures/exch_fluxes_galMut.fig'])
%saveas(fig,['../results/figures/exch_fluxes_galMut.pdf'])
%Visualize results for wild-type model
matrix = abs(exchangesWT);
figure
plot(matrix(:,2),matrix(:,3),matrix(:,2),matrix(:,4),matrix(:,2),matrix(:,5),matrix(:,2),matrix(:,6),matrix(:,2),matrix(:,7),matrix(:,2),matrix(:,8),'LineWidth',4);
xlabel('Oxygen uptake [mmol/gDw h]','FontSize',18)
ylabel('Exchange flux [mmol/gDw h] - [1/h]','FontSize',18)
ylim([0 4])
legend(model.rxnNames(indexes(3:end)))
%title('WT galactose')
%production envelope 
 mutant.lb(indexesM(2)) = -30;
 mutant.ub(indexesM(2)) = 0;
 mutant.lb(indexesM(3)) = maxGM;

 model.lb(indexesM(2)) = -30;
 model.ub(indexesM(2)) = 0;
 model.lb(indexesM(3)) = maxG;
% oxygen = [-10 -2.5 -2 -1.5];
% pos =  find(strcmp(mutant.rxnNames,'galactitol exchange'));
% mutant = setParam(mutant,'obj',pos,1);
% model = setParam(model,'obj','r_4041',1);
% for j=1:4
%     mutant = setParam(mutant,'obj','r_4041',1);
%     model = setParam(model,'obj','r_4041',1);
% 
%      model.lb(indexesM(2)) = oxygen(j);
%      mutant.lb(indexesM(2)) = oxygen(j);
% 
% 
%      solMT = solveLP(mutant);
%      printFluxes(mutant,solMT.x;
%      maxGM = solMT.x(find(mutant.c));
% 
%     mutant = setParam(mutant,'obj',pos,1);
%     model = setParam(model,'obj',pos,1);
% for i=1:iterations+1
%     mutant.lb(indexesM(3)) = (1-((i-1)/iterations))*maxGM;
%     solM = solveLP(mutant);
%     gRateM(i) = solM.x(indexesM(3))/(0.342*solM.x(indexesM(1)));
%     prodM(i) = 0.182*solM.x(pos)/(0.342*solM.x(indexesM(1)));
% end
% plot(-gRateM,-prodM,'LineWidth',4);
% hold on
% end
% 
% mutant.lb(indexesM(2)) = -1000;
% mutant.ub(indexesM(2)) = 0;
% mutant.lb(indexesM(3)) = 0;
% mutant = changeMedia_batch(mutant,'lactose exchange',1);
% mutant = setParam(mutant,'obj','r_4041',1);
% 
% solMG = solveLP(mutant,1);
% 
% oxMut = removeGenes(model,{'xyl1' 'xyl1_2'},false,false,true);
% oxMut = changeMedia_batch(oxMut,'lactose exchange',1);
% oxMut = setParam(oxMut,'obj','r_4041',1);
% oxMut = setParam(oxMut,'lb','r_4041',0);
% solMO = solveLP(oxMut,1);
% 
% diary '../results/CintGEM_lactose_flux_leloir_vs_oxRed_ATP.txt'
% followChanged(model,solMG.x, solMO.x,10,0.05,0.01,'ATP')
% diary off
% 
% diary '../results/CintGEM_lactose_flux_leloir_vs_oxRed_NADH.txt'
% followChanged(model,solMG.x, solMO.x,10,0.05,0.01,'NADH')
% diary off
% 
% diary '../results/CintGEM_lactose_flux_leloir_vs_oxRed_NADPH.txt'
% followChanged(model,solMG.x, solMO.x,10,0.05,0.01,'NADPH')
% diary off

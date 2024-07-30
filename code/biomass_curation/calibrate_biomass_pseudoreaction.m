function modelMod=calibrate_biomass_pseudoreaction(model)
Ptot = 0.438; %average across chemostats in g protein / gCDW

%we have checked that all the metabolites composing the biomass
%reaction are all consistent with those in yeastGEM
comps = getBiomassComponents;
%Get main fractions:
[Pbase,X] = getFraction(model,comps,'P',0);
[Cbase,X] = getFraction(model,comps,'C',X);
[Rbase,X] = getFraction(model,comps,'R',X);
[Dbase,X] = getFraction(model,comps,'D',X);
[Lbase,X] = getFraction(model,comps,'L',X);
%The sum Protein + carbohydrates + lipid backbones + RNA + DNA 
%accounts for 99.75% of the biomass in the model

%Now rescale biomass according to the provided experimental value of total
%protein content. This will rescale carbs and lipids using the same
%proportions as in the base biomass equation
Ctot = Cbase + (Pbase - Ptot)*Cbase/(Lbase+Cbase);
Ltot = Lbase + (Pbase - Ptot)*Lbase/(Lbase+Cbase);
%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fL = Ltot/Lbase;
%Change compositions:
modelMod = rescalePseudoReaction(model,'protein',fP);
modelMod = rescalePseudoReaction(modelMod,'carbohydrate',fC);
%If model contain SLIMER reactions (separate pseudoreactions for
%lipid chains and backbones
modelMod = rescalePseudoReaction(modelMod,'lipid backbone',fL);
modelMod = rescalePseudoReaction(modelMod,'lipid chain',fL);
%Check how stoichiometries have changed for each of the biomass components
% constructEquations(modelMod,posBiomass)
% constructEquations(model,posProt)
% constructEquations(modelMod,posProt)
% constructEquations(model,posCarb)
% constructEquations(modelMod,posCarb)
% constructEquations(model,posLip)
% constructEquations(modelMod,posLip)
%recompute the sum of mass fractions (Protein + carbohydrates + lipid backbones + RNA + DNA) 
[~,X] = getFraction(modelMod,comps,'P',0);
[~,X] = getFraction(modelMod,comps,'C',X);
[~,X] = getFraction(modelMod,comps,'R',X);
[~,X] = getFraction(modelMod,comps,'D',X);
[~,X] = getFraction(modelMod,comps,'L',X);
disp(['Calibrated biomass composition accounts for ' num2str(X) '% of gCDW'])
end
function temp = changePOratio(model,coeff)
temp = model;
rxn(1) = find(strcmp(model.rxns,'r_0438'));%complexIV (half alpha protons are exported)
rxn(2) = find(strcmp(model.rxns,'r_0439'));%complexIII (twice alpha protons are exported)
rxn(3) = find(strcmp(model.rxns,'r_5195'));
met(1) = find(strcmp(model.mets,'s_0794')); %cytoplasmic
met(2) =  find(strcmp(model.mets,'s_0799'));%mito
%Change POratio:
%originalRatio = abs(temp.S(rxn(3),met(2))/temp.S(rxn(2),met(2)));
temp.S(met(2),rxn(2)) = -coeff; 
temp.S(met(1),rxn(2)) = 2*coeff;
temp.S(met(2),rxn(1)) = -coeff; 
temp.S(met(1),rxn(1)) = 0.5*coeff;
temp.S(met(2),rxn(1)) = -coeff; 
temp.S(met(1),rxn(1)) = 0.5*coeff;
%complex I
temp.S(met(2),rxn(3)) = -2.5*coeff;
temp.S(met(1),rxn(3)) = (2.5*coeff)-1;
end
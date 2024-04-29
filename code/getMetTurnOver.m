function [metTmat] = getMetTurnOver(solutions,model)
S = model.S;
S = abs(S);
solu = abs(solutions);
metTmat = 0.5*(S*solu);
end
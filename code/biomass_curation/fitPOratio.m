function POratio = fitPOratio(model)
%Load chemostat data:
fid = fopen('../../data/chemostatData_glucose.txt','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
%exp_data = [exp_data{1} exp_data{2}];
fclose(fid);
%POratios to span:
disp('Estimating POratio:')
%POratio = 5:5:100;
%1st iteration:
%POratio = iteration(model,POratio,exp_data);
%2nd iteration:
%POratio = iteration(model,POratio-10:1:POratio+10,exp_data);
%3rd iteration:
rxn(1) = find(strcmp(model.rxns,'r_0438'));%complexIV (half alpha protons are exported)
rxn(2) = find(strcmp(model.rxns,'r_0439'));%complexIII (twice alpha protons are exported)
rxn(3) = find(strcmp(model.rxns,'r_5195'));
met(1) = find(strcmp(model.mets,'s_0794')); %cytoplasmic
met(2) =  find(strcmp(model.mets,'s_0799'));%mito

alpha  = abs(full(model.S(met(2),rxn(2))));
alphas = (alpha-0.5*alpha):0.01:(alpha+0.5*alpha);
%alphas = 2:0.01:3;
[POratio,error] = iteration(model,alphas,exp_data);
%If verbose output is not required, then the only displayed value is the optimal one
disp(['Fitted POratio = ' num2str(POratio) ' -> Error = ' num2str(error)])

%Plot fit:
mod_data = simulateChemostat(model,exp_data,1,POratio);
mod_data = abs(mod_data);
figure
hold on
cols = [0,1,0;0,0,1;1,0,0];
b    = zeros(1,length(exp_data(1,:))-1);
for i = 1:length(exp_data(1,:))-1
    b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',2);
    plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:))
end
title('POratio fitting for growth on glucose minimal media')
xlabel('Dilution rate [1/h]')
ylabel('Exchange fluxes [mmol/gDWh]')
legend(b,'Glucose consumption','O2 consumption','CO2 production','Location','northwest')
hold off
disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [POratio,error] = iteration(model,POratio,exp_data)

fitting = ones(size(POratio))*1000;

for i = 1:length(POratio)
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model,exp_data,1,POratio(i));
    %mod_data = abs(mod_data(:,[3,4]));
    %exp_data = abs(exp_data(:,[3,4]));
    expData = abs(exp_data);%(:,[3,4]));
    %disp(expData)
    %disp(mod_data)
    R          = (abs(mod_data) - expData)./expData;
    %R          = abs((abs(mod_data(3)) - abs(exp_data(3))))/abs(exp_data(3));
    %fitting(i) = R*100;
    R = R(:,[2,4]);
    fitting(i) = sqrt(sum(sum(R.^2)));
    %disp(['POratio = ' num2str(POratio(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[error,best] = min(fitting);

if best == 1 || best == length(POratio)
    error('POratio found is sub-optimal: please expand POratio search bounds.')
else
    POratio   = POratio(best);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
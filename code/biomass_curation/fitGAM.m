function GAM = fitGAM(model)
% fitGAM
%
% GAM = fitGAM(model,verbose)
% Returns a fitted GAM for the metabolic model (non-GECKO).
%
% Usage: GAM = fitGAM(model)
%

%Load chemostat data:
fid = fopen('../../data/chemostatData_glucose.txt','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
%exp_data = [exp_data{1} exp_data{2}];
fclose(fid);
%GAMs to span:
disp('Estimating GAM:')
GAM = 5:5:100;
%1st iteration:
GAM = iteration(model,GAM,exp_data);
%2nd iteration:
GAM = iteration(model,GAM-10:1:GAM+10,exp_data);
%3rd iteration:
[GAM,error] = iteration(model,GAM-1:0.1:GAM+1,exp_data);
%If verbose output is not required, then the only displayed value is the optimal one
disp(['Fitted GAM = ' num2str(GAM) ' -> Error = ' num2str(error)])

%Plot fit:
mod_data = simulateChemostat(model,exp_data,GAM);
mod_data = abs(mod_data);
figure
hold on
cols = [0,1,0;0,0,1;1,0,0];
b    = zeros(1,length(exp_data(1,:))-1);
for i = 1:length(exp_data(1,:))-1
    b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',2);
    plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:))
end
title('GAM fitting for growth on glucose minimal media')
xlabel('Dilution rate [1/h]')
ylabel('Exchange fluxes [mmol/gDWh]')
legend(b,'Glucose consumption','O2 consumption','CO2 production','Location','northwest')
hold off
disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GAM,error] = iteration(model,GAM,exp_data)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model,exp_data,GAM(i));
    R          = (abs(mod_data) - exp_data)./exp_data;
    %R          = (abs(mod_data(2)) - exp_data(2))/exp_data(2);
    %fitting(i) = R*100;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(exp_data)
    disp(mod_data)
    disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[error,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM   = GAM(best);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

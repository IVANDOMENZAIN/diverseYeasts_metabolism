function NGAM = fitNGAM(model)
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
disp('Estimating NGAM:')
NGAM = 0:0.1:2;
%1st iteration:
NGAM = iteration(model,NGAM,exp_data);
%2nd iteration:
[NGAM,error] = iteration(model,NGAM-0.1:0.01:NGAM+0.1,exp_data);
%If verbose output is not required, then the only displayed value is the optimal one
disp(['Fitted NGAM = ' num2str(NGAM) ' -> Error = ' num2str(error)])

%Plot fit:
mod_data = simulateChemostat(model,exp_data,NGAM);
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

function [NGAM,error] = iteration(model,NGAM,exp_data)

fitting = ones(size(NGAM))*1000;

for i = 1:length(NGAM)
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model,exp_data,[],[],NGAM(i));
    R          = (abs(mod_data) - exp_data)./exp_data;
    %R          = (abs(mod_data(2)) - exp_data(2))/exp_data(2);
    %fitting(i) = R*100;
    R = R(:,[2,4]);
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['NGAM = ' num2str(NGAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[error,best] = min(fitting);

if best == 1 %|| best == length(NGAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    NGAM   = NGAM(best);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

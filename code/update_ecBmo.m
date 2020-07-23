current = pwd;
%Clone GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('pull')
%Locate the correct branch
%git('stash')
git('checkout fix/updateDatabases')
cd ..
model = importModel('../models/blastobotrys_mokoenaii.xml');
cd GECKO/geckomat/get_enzyme_data
[uni,kegg] = updateDatabases;
cd ../../..
toolbox = 'COBRA';
name='ecBmo';
version='1.0';
pause
%Replace scripts in GECKO:
scripts = dir('specific_scripts');
for i = 1:length(scripts)
    script = scripts(i).name;
    if contains(script,'.m')
        fullName   = ['specific_scripts/' script];
        %Retrieve script path within GECKO
        GECKO_path = dir(['GECKO/' script]);
        if ~isempty(GECKO_path)
            GECKO_path = GECKO_path.folder;
            %Replace script in GECKO in its container subfolder
            copyfile(fullName,GECKO_path)
        end
    end
end

[ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version);
%COpy outputs to repo
cd ../..
copyfile(['GECKO/models/' name],['../models/' name])
save(['../models/' name '/ecModel.mat'],'ecModel')
save(['../models/' name '/ecModel_batch.mat'],'ecModel_batch')
%xylose growth
%
% Ivan Domenzain. 2020-05-06
%
current = pwd;
%Clone GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('pull')
%Locate the correct branch
%git('stash')
git('checkout feat-add_utilities')
clc
cd ..
%Retrieve model names
originFolder = '../Reconstruction_script/ModelFiles/xml';
fileNames   = dir(originFolder);
namesT      = [];
xylRxnT     = [];
araRxnT     = [];
gRates_glcT = [];
gRates_xylT = [];
gRates_araT = [];

for i=1:length(fileNames)
    cd(current)
    file = fileNames(i).name;
    if contains(file,'.xml')
        modelName = file(1:end-4);
        if startsWith(modelName,'yH')
            pos = strfind(modelName,'_');
            modelName = modelName((pos(1)+1):end);
        end    
        namesT  = [namesT;{modelName}];
        gRate_xyl = 0;
        gRate_glc = 0;
        gRate_ara = 0;
        disp(['******* ' modelName ' *******'])
        model = importModel([originFolder '/' file]);
        %Set growth as objective function
        model = setParam(model,'obj',{'r_2111'},1);
        gRate_glc = batchGrowth(model,'r_1714');
        disp(['The gRate on D-glucose is: ' num2str(gRate_glc)])
        %Growth on Xylose
        xylRxn = find(strcmpi(model.rxns,'r_1718'));
        if ~isempty(xylRxn)
            xylRxnT   = [xylRxnT; xylRxn];
            gRate_xyl = batchGrowth(model,'r_1718');
            disp(['The gRate on D-xylose is: ' num2str(gRate_xyl)])
        else
            xylRxnT = [xylRxnT; 0];
        end
        %Growth on arabinose
        araRxn = find(strcmpi(model.rxns,'r_1878'));
        if ~isempty(araRxn)
            araRxnT   = [araRxnT; araRxn];
            gRate_ara = batchGrowth(model,'r_1878');
            disp(['The gRate on L-arabinose is: ' num2str(gRate_ara)])
        else
            araRxnT = [araRxnT; 0];
        end
        gRates_xylT = [gRates_xylT;gRate_xyl];
        gRates_glcT = [gRates_glcT;gRate_glc];
        gRates_araT = [gRates_araT;gRate_ara];
    end
    disp(' ')
end
T = table(namesT,xylRxnT,gRates_glcT,gRates_xylT,gRates_araT);
writetable(T,'results_xylose.txt','Delimiter','\t','QuoteStrings',false)

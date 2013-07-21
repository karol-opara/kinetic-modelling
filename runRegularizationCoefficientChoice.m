function runRegularizationCoefficientChoice(name)

load saveExperimentalData20130319_8of13experiments_NRTLvalidation
poorData = data{4};
goodData = data{2};
%plotConcentrations(poorData.timeX, poorData.xml,poorData.timeY, poorData.yml,poorData.timeZ, poorData.zml, '-');
%plotConcentrations(goodData.timeX, goodData.xml,goodData.timeY,goodData.yml,goodData.timeZ, goodData.zml, '-');

savefilename = ['Results/' 'save_RegularizationCoefficients_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];

lambdaZero = [1, 0];

type = 'batch';

tic

p = {'log', 0.5,  1,  2}; % inner
q = {NaN, 'log', 0.5,  1,  2}; % outer loss (NaN means using only inner loss)

expErrCoef=1.05;

pm = cell(length(p),length(q));
errMin = cell(length(p),length(q));
errMinReg = cell(length(p),length(q));
lambda = cell(length(p),length(q));
for i=1:length(p)
    parfor j = 1:length(q)
        pm{i,j} = EstimateKineticModel(poorData,p{i},q{j},type,lambdaZero);
        errMin{i,j} = pm{i,j}.optimizerOutput.solutions.bestever.f;
        
        expectedError = errMin{i,j}*expErrCoef;
        options = optimset('fmincon');
        options.TolFun = 0.1;
        [lambda{i,j},errMinReg{i,j}] = fmincon(@EstimateRegularizedModel,0,[],[],[],[],[0],[],[],options,...
            poorData,p{i},q{j},type, expectedError);
        
        relErr = abs(errMinReg{i,j}/errMin{i,j}-expectedError);
        if(relErr>0.05)
            warning('runRegularizationCoefficientChoice:RelErr',['too large relative error ' num2str(relErr)]);
        end
    end
    fprintf('.');
end



save(savefilename);
disp(['Saved as ' savefilename]);
toc
end

function err = EstimateRegularizedModel(lambdaReg, data,p,q,type, expectedError)
model = EstimateKineticModel(data,p,q,type,[1, lambdaReg]);
error = model.optimizerOutput.solutions.bestever.f;
err = abs(error-expectedError);
end
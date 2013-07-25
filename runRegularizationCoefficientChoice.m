function runRegularizationCoefficientChoice(name)
loadData = true;

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
relativeAccuracy = 0.005;

pm = cell(length(p),length(q));
errMin = zeros(length(p),length(q));
errMinReg = zeros(length(p),length(q));
lambda = zeros(length(p),length(q));

if (loadData)
    load('Results/save_RegularizationCoefficients_2013-07-23_092626_1e5Dim_LBFGS_PoorData_NonregularizedErrors');
else
    options.MaxFunEvals = 6*1e5;
    for i=1:length(p)
        parfor j = 1:length(q)
            pm{i,j} = EstimateKineticModel(poorData,p{i},q{j},type,lambdaZero,options);
            errMin(i,j) = pm{i,j}.optimizerOutput.solutions.bestever.f;
        end
        fprintf('.');
    end
    save([savefilename '_NonregularizedErrors']);
    disp(['Saved nonregularized errors as ' savefilename  '_NonregularizedErrors']);
    fprintf('\n');
end

opts = cmaes('defaults');
opts.MaxFunEvals = 1e3;
opts.MaxIter = Inf;
opts.LBounds = 0;
opts.UBounds = 100; % i.e. over 85 = 0.05 [accuracy] * max(max(errMin)) [no regularization erro] / 0.18 [smallest norm of reaction rates in literature]
opts.SaveVariables = 'off';
opts.DispFinal = 0;
opts.DispModulo = 0;
opts.SaveFilename = 0;
opts.LogFilenamePrefix = 0;
opts.LogTime = 1;
opts.LogModulo = Inf;
warning('off','cmaes:logging');
opts.Restarts = 0;
cmaesOut = [];

for i=1:length(p)
    parfor j = 1:length(q)
        expectedError = errMin(i,j)*expErrCoef;
        jOpts = opts;
        jOpts.StopFitness = relativeAccuracy*errMin(i,j);
        %[lambda(i,j),errMinReg(i,j)] = fminlbfgs(@EstimateRegularizedModel,rand(),options,poorData,p{i},q{j},type, expectedError);
        [lambda(i,j), errMinReg(i,j), ~, flag, cmaesOut(i,j), ~]= cmaes('prvEstimateRegularizedModel', rand() * (opts.LBounds+opts.UBounds)/4, [], jOpts, poorData, p{i},q{j},type, expectedError);
        if (any(strcmp(flag,'fitness')) == 0)
            warning('runRegularizationCoefficientChoice:CmaesStopFlag',['Unexpected stop flag: ' flag]);            
        end
        relErr = abs(errMinReg(i,j)/errMin(i,j)-expectedError);
        if(relErr>0.05)
            warning('runRegularizationCoefficientChoice:RelErr',['too large relative error ' num2str(relErr)]);
        end
    end
    save([savefilename '_partial']);
    fprintf('.');
end

save(savefilename);
disp(['Saved as ' savefilename]);
toc
end



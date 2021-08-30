function runRegularizationCoefficientChoice(name)
loadData = false;
optimization = false;

load saveExperimentalData20130319_8of13experiments_NRTLvalidation
poorData = data{4};
goodData = data{2};
%plotConcentrations(poorData.timeX, poorData.xml,poorData.timeY, poorData.yml,poorData.timeZ, poorData.zml, '-');
%plotConcentrations(goodData.timeX, goodData.xml,goodData.timeY,goodData.yml,goodData.timeZ, goodData.zml, '-');

savefilename = strjoin(['Results/' 'save_RegularizationCoefficients_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name], '');
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

if (optimization)
    for rep = 1:12
        if (loadData)
            load('Results/save_RegularizationCoefficients_2013-07-23_092626_1e5Dim_LBFGS_PoorData_NonregularizedErrors');
            warning('runRegularizationCoefficient','Data on the unregularized errors loaded instead of computing');
        else
            % options.MaxFunEvals = 6*1e4;
            options = struct();
            % Instead of writing two nested loops for i=1:length(p)
            % for j = 1:length(q)} we write a single loop, which
            % parallelizes better with parfor
            parfor k = 1:(length(p) * length(q))
                i = 1+floor((k-1)/length(q));
                j = 1+mod(k-1, length(q));
                pmk{k} = EstimateKineticModel(poorData,p{i},q{j},type,lambdaZero,options,'madDE');
                errMink(k) = pmk{k}.optimizerOutput.best_val;
            end
            for k = 1:(length(p) * length(q))
                i = 1+floor((k-1)/length(q));
                j = 1+mod(k-1, length(q));
                pm{i,j} = pmk{k};
                errMin(i,j) = errMink(k);
            end
            
            %     for i=1:length(p)
            %         parfor j = 1:length(q)
            %             pm{i,j} = EstimateKineticModel(poorData,p{i},q{j},type,lambdaZero,options,'madDE');
            %             errMin(i,j) = pm{i,j}.optimizerOutput.best_val;
            %         end
            %     end
            save(strjoin([savefilename '_NonregularizedErrors_run_' num2str(rep)], ''));
            disp(strjoin(['Saved nonregularized errors as ' savefilename  '_NonregularizedErrors_run_' num2str(rep)],''));
            fprintf('\n');
        end
    end
end


fprintf('Iterating over relative lambdas\n');
relativeLambdas = [1 0.5 0.2 0.1 0.05 0.02 0.01]; % [0.7 0.07] [2 5] [7 10]
for i = 1:length(relativeLambdas)
    lambda = [1 relativeLambdas(i)];
    runBenchmarkingExperiment('RegularizationCoefficientChoice', ...
        0.05, [1 0 0 0 0 0], 0.01, 'Oh', lambda, 12, NaN, false);
    save(strjoin([savefilename '_partial'], ''));
    fprintf('.');
end


save(savefilename);
disp(strjoin(['Saved as ' savefilename],''));
toc
end

function dispProgress(m)
fprintf([num2str(m) ' ']);
end


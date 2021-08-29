function processRegularizationCoefficientChoice(name, wd)
if (nargin < 1)
    error('Experiment name cannot be empty')
end
threshold = 1;
[err, relErr, lambda2] = gerErrors(threshold, wd);

bestErr = Inf(4,5);
bestLambda = zeros(4,5);
for i=1:size(relErr,1)
    for j=1:size(relErr,2)
        for l = 1:size(relErr,3)
            if(relErr(i,j,l)<bestErr(i,j))
                bestErr(i,j) = relErr(i,j,l);
                bestLambda(i,j) = lambda2(l);
            end
        end
    end
end
bestErrPercentage = round(100*bestErr)
bestLambda

savefilename = ['Results/' 'save_RegularizationCoefficientChoice_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];
save(savefilename);
end

function [err, relErr, lambda2] = gerErrors(relErrTreshold, wd)
err=[];

fns = {'save_2021-08-28_120803BenchmarkingExperiment_12_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice_minErr_001_partial.mat',...
'save_2021-08-28_133500BenchmarkingExperiment_12_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-28_143409BenchmarkingExperiment_12_Repetetive_Fits_105_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-28_180823BenchmarkingExperiment_12_Repetetive_Fits_102_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-28_190817BenchmarkingExperiment_12_Repetetive_Fits_101_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-28_205151BenchmarkingExperiment_12_Repetetive_Fits_1005_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-28_215318BenchmarkingExperiment_12_Repetetive_Fits_1002_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_084714BenchmarkingExperiment_12_Repetetive_Fits_1001_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_094934BenchmarkingExperiment_12_Repetetive_Fits_10005_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_131253BenchmarkingExperiment_12_Repetetive_Fits_12_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_141410BenchmarkingExperiment_12_Repetetive_Fits_15_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_194844BenchmarkingExperiment_12_Repetetive_Fits_17_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
'save_2021-08-29_204838BenchmarkingExperiment_12_Repetetive_Fits_110_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat'};

for i = 1:length(fns)

[err(:,:,i), relErr(:,:,i), lambda2(i)] = getError(fullfile(wd, 'Results', fns{i}),...
    relErrTreshold);

end


% [err(:,:,1), relErr(:,:,1), lambda2(1)] = getError('Results\save_2013-07-29_140735BenchmarkingExperiment_14_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,2), relErr(:,:,2), lambda2(2)] = getError('Results\save_2013-07-30_010755BenchmarkingExperiment_14_Repetetive_Fits_105_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,3), relErr(:,:,3), lambda2(3)] = getError('Results\save_2013-07-30_123147BenchmarkingExperiment_14_Repetetive_Fits_102_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,4), relErr(:,:,4), lambda2(4)] = getError('Results\save_2013-07-31_001956BenchmarkingExperiment_14_Repetetive_Fits_101_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,5), relErr(:,:,5), lambda2(5)] = getError('Results\save_2013-07-31_122905BenchmarkingExperiment_14_Repetetive_Fits_1005_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,6), relErr(:,:,6), lambda2(6)] = getError('Results\save_2013-08-01_011632BenchmarkingExperiment_14_Repetetive_Fits_1002_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
% [err(:,:,7), relErr(:,:,7), lambda2(7)] = getError('Results\save_2013-08-01_144049BenchmarkingExperiment_14_Repetetive_Fits_1001_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);

end

function [err, relErr, lambda2] = getError(filename, relErrTreshold)
load(filename);
lambda2=lambda(2);
[ii,jj,ll] = size(models);
err = zeros(jj,ll);
relErr=zeros(jj,ll);

kRefOh = [0.0500
    0.1100
    0.2150
    1.2280
    0.2420
    0.0070];

for j=1:jj
    for l=1:ll
        relError=NaN(1,ii);
        for i=1:ii
            relError(i) = getRelativeError(kRefOh, models{i,j,l}.k);
            if (relError(i) < relErrTreshold)
                err(j,l) = err(j,l) + 1;
            end
        end
        relError = sort(relError);
        relErr(j,l) = mean(relError); % relError(1:10)
    end
end
end

function relErr = getRelativeError(kOriginal,kEstimated)
relErr = max(abs(kEstimated-kOriginal)./kOriginal);
end
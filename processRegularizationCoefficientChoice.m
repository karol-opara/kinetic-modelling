function processRegularizationCoefficientChoice(name)
if (nargin < 1)
    error('Experiment name cannot be empty')
end
threshold = 1;
[err, relErr, lambda2] = gerErrors(threshold);

bestErr = Inf(4,5);
bestLambda = zeros(4,5);
for i=1:4
    for j=1:5
        for l = 1:7
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

function [err, relErr, lambda2] = gerErrors(relErrTreshold)
err=[];
[err(:,:,1), relErr(:,:,1), lambda2(1)] = getError('Results\save_2013-07-29_140735BenchmarkingExperiment_14_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,2), relErr(:,:,2), lambda2(2)] = getError('Results\save_2013-07-30_010755BenchmarkingExperiment_14_Repetetive_Fits_105_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,3), relErr(:,:,3), lambda2(3)] = getError('Results\save_2013-07-30_123147BenchmarkingExperiment_14_Repetetive_Fits_102_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,4), relErr(:,:,4), lambda2(4)] = getError('Results\save_2013-07-31_001956BenchmarkingExperiment_14_Repetetive_Fits_101_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,5), relErr(:,:,5), lambda2(5)] = getError('Results\save_2013-07-31_122905BenchmarkingExperiment_14_Repetetive_Fits_1005_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,6), relErr(:,:,6), lambda2(6)] = getError('Results\save_2013-08-01_011632BenchmarkingExperiment_14_Repetetive_Fits_1002_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
[err(:,:,7), relErr(:,:,7), lambda2(7)] = getError('Results\save_2013-08-01_144049BenchmarkingExperiment_14_Repetetive_Fits_1001_RelativeLambda_RegularizationCoefficientChoice',relErrTreshold);
end

function [err, relErr, lambda2] = getError(filename, relErrTreshold)
load(filename);
lambda2=lambda(2);
[ii,jj,ll] = size(models);
err = zeros(jj,ll);
relErr=zeros(jj,ll);

for j=1:jj
    for l=1:ll
        relError=NaN(1,ii);
        for i=1:ii
            relError(i) = getRelativeError(models{i,j,l}.data.k, models{i,j,l}.k);
            if (relError(i) < relErrTreshold)
                err(j,l) = err(j,l) + 1;
            end
        end
        relError = sort(relError);
        relErr(j,l) = mean(relError(1:10));
    end
end
end

function relErr = getRelativeError(kOriginal,kEstimated)
relErr = max(abs(kEstimated-kOriginal)./kOriginal);
end
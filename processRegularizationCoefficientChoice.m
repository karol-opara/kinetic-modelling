function processRegularizationCoefficientChoice
[err, lambda2] = getError('Results\save_2013-07-28_001607BenchmarkingExperiment_14_Repetetive_Fits_ 01_RelativeLambda_RegularizationCoefficientChoice')

[err, lambda2] = getError('Results\save_2013-07-28_121843BenchmarkingExperiment_14_Repetetive_Fits_005_RelativeLambda_RegularizationCoefficientChoice')
end

function [err, lambda2] = getError(filename)
relErrTreshold = 0.050;

load(filename);
lambda2=lambda(2);
[ii,jj,ll] = size(models);
err = zeros(jj,ll);
for i=1:ii
    for j=1:jj
        for l=1:ll
            relErr = getRelativeError(models{i,j,l}.data.k, models{i,j,l}.k);
            if (relErr < relErrTreshold)
                err(j,l) = err(j,l) + 1/ii;
            end
        end
    end
end
err= round(err*100);
end

function relErr = getRelativeError(kOriginal,kEstimated)
relErr = max(abs(kEstimated-kOriginal)./kOriginal);
end
function processRegularizationCoefficientChoice(name, wd)
if (nargin < 1)
    error('Experiment name cannot be empty')
end
threshold = 1;
[err, relErr, lambda2, errL2, stdErrL2] = gerErrors(threshold, wd);

figure()
bestErr = Inf(4,5);
bestLambda = zeros(4,5);
k = 1;
for i=1:size(errL2,1)
    for j=1:size(errL2,2)
        subplot(size(errL2,1), size(errL2, 2), k);
        k = k + 1;
        e2 = NaN(size(errL2,3),1);
        s2 = NaN(size(errL2,3),1);
        for l = 1:size(errL2,3)
            e2(l) = errL2(i,j,l);
            s2(l) = stdErrL2(i,j,l);
            
            if(errL2(i,j,l)<bestErr(i,j))
                bestErr(i,j) = errL2(i,j,l);
                bestLambda(i,j) = lambda2(l);
            end
        end
        [~, is] = sort(lambda2);
        h = errorbar(lambda2(is), e2(is), s2(is));
        xlim([-0.1 0.6]);
        
        y = ylim; % current y-axis limits
        hl = line([bestLambda(i,j) bestLambda(i,j)], [y(1) y(2)]);
        set(hl, 'LineStyle', ':');
        ylim(y);
    end
end
bestErrPercentage = round(100*bestErr);
bestLambda

savefilename = ['../results/' 'save_RegularizationCoefficientChoice_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];
save(savefilename);
disp(['Saved results as ' savefilename]);
end

function [err, relErr, lambda2, errL2, stdErrL2] = gerErrors(relErrTreshold, wd)
err=[];

% We don't want too large regularization coefficients.
%'save_2021-09-02_124353BenchmarkingExperiment_12_Repetetive_Fits_12_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
%'save_2021-08-30_150609BenchmarkingExperiment_12_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
%'save_2021-08-31_081100BenchmarkingExperiment_12_Repetetive_Fits_107_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...

fns = {'save_2021-08-30_160310BenchmarkingExperiment_12_Repetetive_Fits_105_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
    'save_2021-08-30_171901BenchmarkingExperiment_12_Repetetive_Fits_102_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
    'save_2021-08-30_181753BenchmarkingExperiment_12_Repetetive_Fits_101_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
    'save_2021-08-30_193602BenchmarkingExperiment_12_Repetetive_Fits_1005_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
    'save_2021-08-30_203627BenchmarkingExperiment_12_Repetetive_Fits_1002_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...save_2021-08-30_215103BenchmarkingExperiment_12_Repetetive_Fits_1001_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
    'save_2021-08-31_090733BenchmarkingExperiment_12_Repetetive_Fits_1007_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat'};

% fns = {'save_2021-08-28_120803BenchmarkingExperiment_12_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice_minErr_001_partial.mat',...
% 'save_2021-08-28_133500BenchmarkingExperiment_12_Repetetive_Fits_11_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-28_143409BenchmarkingExperiment_12_Repetetive_Fits_105_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-28_180823BenchmarkingExperiment_12_Repetetive_Fits_102_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-28_190817BenchmarkingExperiment_12_Repetetive_Fits_101_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-28_205151BenchmarkingExperiment_12_Repetetive_Fits_1005_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-28_215318BenchmarkingExperiment_12_Repetetive_Fits_1002_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_084714BenchmarkingExperiment_12_Repetetive_Fits_1001_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_094934BenchmarkingExperiment_12_Repetetive_Fits_10005_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_131253BenchmarkingExperiment_12_Repetetive_Fits_12_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_141410BenchmarkingExperiment_12_Repetetive_Fits_15_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_194844BenchmarkingExperiment_12_Repetetive_Fits_17_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat',...
% 'save_2021-08-29_204838BenchmarkingExperiment_12_Repetetive_Fits_110_RelativeLambda_RegularizationCoefficientChoice_minErr_001.mat'};

for i = 1:length(fns)
    
    [err(:,:,i), relErr(:,:,i), lambda2(i), errL2(:,:,i), stdErrL2(:,:,i)] = getError(fullfile(wd, '../data/regularization', fns{i}),...
        relErrTreshold);
    
end
end

function [err, relErr, lambda2, errL2, stdErrL2] = getError(filename, relErrTreshold)
load(filename, 'models', 'lambda');
lambda2=lambda(2);
[ii,jj,ll] = size(models);
err = zeros(jj,ll);
errL2 = zeros(jj, ll);
stdErrL2 = zeros(jj, ll);
relErr=zeros(jj,ll);

kRefOh = [0.0500
    0.1100
    0.2150
    1.2280
    0.2420
    0.0070];

for j=1:jj
    for l=1:ll
        errorL2 = NaN(1,ii);
        relError=NaN(1,ii);
        for i=1:ii
            errorL2(i) = sqrt(sum((kRefOh - models{i,j,l}.k).^2));
            relError(i) = getRelativeError(kRefOh, models{i,j,l}.k);
            if (relError(i) < relErrTreshold)
                err(j,l) = err(j,l) + 1;
            end
        end
        errL2(j,l) = mean(errorL2);
        stdErrL2(j,l) = std(errorL2);
        relError = sort(relError);
        relErr(j,l) = mean(relError); % relError(1:10)
    end
end
end

function relErr = getRelativeError(kOriginal,kEstimated)
relErr = max(abs(kEstimated-kOriginal)./kOriginal);
end
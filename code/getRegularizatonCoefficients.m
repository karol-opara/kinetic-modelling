function [lambdas, lossMult] = getRegularizatonCoefficients(wd, i, j)
if (nargin == 3)
    lambdas = NaN;
    lossMult = lossFunctionOrderMultiplier(i,j);
elseif (nargin == 1)
    lossMult = NaN;
    %load([wd '/../data/save_RegularizationCoefficientChoice_2013-08-06_162711_14runsOhsData.mat'])
    load([wd '/../data/regularization/save_RegularizationCoefficientChoice_2021-09-06_184402_relativeLambdaM.mat'], 'bestLambda')
    lambdas = bestLambda;
end
end

function lossMult = lossFunctionOrderMultiplier(i,j)
errMins = NaN(4,5);

% Load the files and compute average min errror
for ii = 1:12
    fn = ['../data/regularization/',...
        'save_RegularizationCoefficients_2021-08-29_235841_relativeLambdaMultipleRuns_NonregularizedErrors_run_', ...
        num2str(ii), '.mat'];
    load(fn, 'errMin')
    %errMin = loadErrMin(fn);
    if ii == 1
        errMins = errMin;
    else
        errMins = (ii-1)/ii * errMins + 1/ii * errMin; % computing arithmetic mean for streaming data
    end
end

lossMult = errMins(i,j);
end
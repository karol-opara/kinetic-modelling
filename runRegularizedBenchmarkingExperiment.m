function runRegularizedBenchmarkingExperiment(name, ... % experiment name (used for identification in output files only)
    dataset,  ... % one of 'Oh' 'Jansri' 'Noureddini' 'Klofutar'
    useLambdas,  ... % should regularization be used?
    useSystematicErr,  ... % should systematic errors be simulated?
    useMinErr, ... % should minimal error be simulated?
    N) % number of runs
if (nargin ~= 6)
    error('Incorrect number of arguments given')
end

lambdas = getRegularizatonCoefficients();
if (useLambdas)
    name = [name '_regularized'];
else
    lambdas = lambdas*0;
    name = [name '_nonregularized'];
end

randomErr = 0.05;
systematicErr = [0 0 0 0 0 0];

if (useSystematicErr)
    systematicErr = [1 0 0 0 0 0];
end

minErr = 0;
if (useMinErr)
    minErr = 0.01;
end

runBenchmarkingExperiment(name, randomErr, systematicErr, minErr, dataset, NaN, N, lambdas);
end

function lambdas = getRegularizatonCoefficients()
load('Results\save_RegularizationCoefficientChoice_2013-08-06_162711_14runsOhsData');
lambdas = bestLambda;
end


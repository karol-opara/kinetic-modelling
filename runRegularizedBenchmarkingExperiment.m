function runRegularizedBenchmarkingExperiment(name, dataset, useLambdas, useSystematicErr, useMinErr, N)

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


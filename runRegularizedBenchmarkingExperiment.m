function runRegularizedBenchmarkingExperiment(name,dataset, useLambdas)
N = 100;
lambdas = getRegularizatonCoefficients();
if (nargin > 2)
    lambdas = lambdas*0;
    name = [name '_nonregularized'];
else
    name = [name '_regularized'];
end    
runBenchmarkingExperiment(name, 0.05, [0 0 0 0 0 0], dataset, NaN, N, lambdas);
end

function lambdas = getRegularizatonCoefficients()
load('Results\save_RegularizationCoefficientChoice_2013-08-06_162711_14runsOhsData');
lambdas = bestLambda;
end


function err = prvEstimateRegularizedModel(lambdaReg, data,p,q,type, expectedError)
model = EstimateKineticModel(data,p,q,type,[1, lambdaReg]);
error = model.optimizerOutput.solutions.bestever.f;
err = (error-expectedError).^2;
end


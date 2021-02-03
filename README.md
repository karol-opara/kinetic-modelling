# kinetic-modelling
Code related to the paper on regularization and concave loss functions for estimation of chemical kinetic models

Integration of a single kinetic model for given initial concentrations and rate constants is implemented in function `SolveKineticModel`. This is the objective function for the optimization routine from `EstimateKineticModel`. 

The benchmarking study is started by calling the `runRegularizedBenchmarkingExperiment` function. Its arguments are the test name (used in filename for storing the output data), the data set, which is the basis of a benchmark ('Oh', 'Jansri', 'Noureddini' or 'Klofutar'), three boolean flags showing whether to use regularization and whether the proportional error should be summed with systematic error and with a minimal error. The number of repeats is the final argument.

A call example:
`runRegularizedBenchmarkingExperiment('test_name', 'Oh', true, true, true, 15)`

Paths to the result files are manually provided in the `processBenchmarkingExperiment` script and should be amended to process new results. Function `processBenchmarkingExperiment` performs a complete analysis of the benchmarking experiment. Function `processLandscapes` generates the error landscapes and plots best-performing models.
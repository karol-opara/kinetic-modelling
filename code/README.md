# kinetic-modelling
Code related to the paper by Karol Opara and Pin Pin Oh in Applied Soft Computing titled "Regularization and concave loss functions for estimation of chemical kinetic models".

Integration of a single kinetic model for given initial concentrations and rate constants is implemented in function `SolveKineticModel`. This is the objective function for the optimization routine from `EstimateKineticModel`. 

The benchmarking study is started by calling the `runRegularizedBenchmarkingExperiment(name, dataset, useLambdas, useSystematicErr, useMinErr, N, compareOptimizers)` function. Its arguments are the test name (used in filename for storing the output data), the data set, which is the basis of a benchmark ('Oh', 'Jansri', 'Noureddini' or 'Klofutar'), three boolean flags showing whether to use regularization and whether the proportional error should be summed with systematic error and with a minimal error. The number of problem instances to be generated (and optimization routine runs). Finally, there is a boolean flag indicating which experiment should be run. Setting it to true compares different optimizers for three selected loss functions as in Table 7 in the paper. Setting the flag to false compares different combinations of the loss function and regularization as in Table 5 in the paper.

After finishing (usually time-consuming) computations, their results are stored as files. This separates the generation of results from their analysis. Paths to the result files are manually provided in the `processBenchmarkingExperiment` script and should be modified if one desires to process newly generated results. Function `processBenchmarkingExperiment` performs a complete analysis of the benchmarking experiment. In particular, it produces Latex tables for presenting the results. 

A comparison of optimizers for three selected loss functions, provided in Table 7 in the paper, was obtained from function `processOptimizationComparison`. The main method reads a cell array of paths to the files with simulation results. Editing of these paths allows for the application of this function for new data.

A call example for loss function and regularization comparison (reproducing Table 5):

`runRegularizedBenchmarkingExperiment('test_name', 'Noureddini', true, true, true, 100, false) ` 

`processBenchmarkingExperiment % Processing of the results (caution: this requires modification of the file paths in the main method to process new results)`

A call example for optimizer comparison (reproducing Table 7):

`runRegularizedBenchmarkingExperiment('test_name', 'Noureddini', true, true, true, 30, true) % Processing of the results (caution: this requires modification of the filepaths in the main method to process new results)` 

`processOptimizationComparison` 


The results of simulations that are discussed in the paper are provided as .mat files so that processing can be started from them. 

The residual landscape plots and trajectories of the best-performing models, presented as Figures 6, 7, 8 and 9 in the paper, are generated by function `processLandscapes`.
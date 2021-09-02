function runBenchmarkingExperiment(name, rndErr, sysErr, minErr, dataset, lambda, N, lambdas, compareOptimizers)

if (nargin < 2)
    rndErr = 0.15;
    sysErr = [1 0 0 0 0 0];
end
if (nargin < 5)
    dataset = 'Oh';
    warning('runBenchmarkingExperiment:NoExperimentID','No experiment ID provided');
end

% if (areTheseToolboxesInstalled({'MATLAB','Parallel Computing Toolbox'}))
%     matlabpool(feature('numCores')-1);
% end

% RunUniqunessExperimentImportanceSampling(randErr, sysErr)
if compareOptimizers
    RunOptimizerComparison(name, rndErr, sysErr, minErr, dataset, N, lambdas);
else
    RunUniqunessExperimentRepetetiveFits(name, rndErr, sysErr, minErr, dataset,lambda,N,lambdas);
end


end


function RunUniqunessExperimentRepetetiveFits(name, rndErr, sysErr, minErr, id,lambda, N,lambdas)
savefilename = ['Results/' 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    'BenchmarkingExperiment_' num2str(N) '_Repetetive_Fits_' num2str(lambda)...
    '_RelativeLambda_' name '_minErr_' num2str(minErr)];
savefilename(ismember(savefilename,' ,.:;!'))=[];

disp(['Repetitive fits']);

%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only 10 repeats');
pnorms = {'log', 0.5, 1, 2};
qnorms = {NaN,'log', 0.5, 1, 2};
%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only NaN norms tried');
plen = length(pnorms);
qlen = length(qnorms);
models = cell(N,plen,plen);
for i = 1:N
    dataN(i) = CreateBenchmarkProblem(rndErr, sysErr, minErr, id);
end
data = dataN(1);

% For better parallelization we flatten three nested loops into a single
% one by (1) obtaining all the indices, (2) running computations in a
% flattened loop, (3) rewriting results to the original data structure
k = 1;
ijrep = [];
for i = 1:plen
    for j = 1:qlen
        for rep = 1:N
            ijrep = [ijrep; i j rep k];
            k = k+1;
        end
    end
end

mk = cell(plen * qlen *N, 1);
k = 1;
parfor k = 1:(plen * qlen *N)
    i = ijrep(k, 1);
    j = ijrep(k, 2);
    rep = ijrep(k, 3);
    
    if (any(any(isnan(lambdas))))
        lambda_reg = [1, lambda*lossFunctionOrderMultiplier(i,j)];
    else
        lambda_reg = [1, lambdas(i,j)*lossFunctionOrderMultiplier(i,j)];
    end
    data = dataN(rep);
    p = pnorms{i};
    q = qnorms{j};
    
    mk{k} = EstimateKineticModel(data,p,q,'batch',lambda_reg,...
        struct(), 'madDE');
    
    fprintf('.');
    % save([savefilename '_partial']);
end

k = 1;
ijrep = [];
for i = 1:plen
    for j = 1:qlen
        for rep = 1:N
            models{rep,i,j} = mk{k};
            k = k+1;
        end
    end
end

save(savefilename);
disp(['Saved as ' savefilename])
end

function lossMult = lossFunctionOrderMultiplier(i,j)
%load('Results/save_RegularizationCoefficients_2013-07-23_092626_1e5Dim_LBFGS_PoorData_NonregularizedErrors');
%warning('Loading old multipliers');

errMins = NaN(4,5);

% Load the files and compute average min errror
for ii = 1:12
    fn = ['Regularization/',...
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


function data = CreateBenchmarkProblem(rndErr, sysErr, minErr, id)
% first
if (strcmp(id,'Oh'))
    k = [0.0500
        0.1100
        0.2150
        1.2280
        0.2420
        0.0070];
elseif (strcmp(id,'Jansri'))
    k = [2.6
        0.248
        1.186
        0.227
        2.303
        0.022];
elseif (strcmp(id,'Noureddini'))
    k = [0.0500
        0.1100
        0.2150
        1.2280
        0.2420
        0.0070];
elseif (strcmp(id,'Klofutar'))
    k = [0.0443
        0.2334
        0.0645
        0.0699
        0.2681
        0.0047];
else
    error('runBenchmarkingExperiment:WrongID','Wrong ID of the reaction rates');
end

conc = [6 0 1 0 0 0];

% totalMolePerLitre = 10; % in simulations and initial submission
totalMolePerLitre = 0.789359594 + 0.057795443 + 0.003075868 + 0.003075868 + 5.002450688; % in revision

z0ml = totalMolePerLitre * conc/sum(conc);
if (strcmp(id,'Oh'))
    timeZ = [2     4     6     8    10    12    20    30    40    50    60];
elseif(strcmp(id,'Jansri'))
    timeZ = [0.5, 1, 3, 5, 7, 9, 12, 15, 18, 20 ];
elseif (strcmp(id,'Noureddini'))
    timeZ = [1 2 3 4 5 6 8 10 15 20 30 50 60 90]; % Nouredini and Zhou (1997)
elseif (strcmp(id, 'Klofutar'))
    timeZ = [2 8 12 16 20 30 50 70 90]; % Pixel-wise reverse engineering from Fig. 4B
end
[t, z] = SolveKineticModel(z0ml,k);
zKinetic = interp1(t,z,timeZ);

%rng('shuffle');

rnd = rndErr* sqrt(pi/2) * randn(length(timeZ),length(z0ml)) + minErr * randn(length(timeZ),length(z0ml));
zml = zKinetic+zKinetic.*rnd+repmat(sysErr,length(timeZ),1);

data = struct('k',k,'z0ml',z0ml,'timeZ',timeZ,'zml',zml,'zKinetic',zKinetic);

end

function RunOptimizerComparison(name, rndErr, sysErr, minErr, id, N,lambdas)
savefilename = ['Results/' 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    '_OptimizationComparison_' id, '_', num2str(N) '_RepetetiveFits_' ...
    name '_minErr_' num2str(minErr)];
savefilename(ismember(savefilename,' ,.:;!'))=[];

disp(['Comparing optimizers']);

%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only 10 repeats');
pnorms = {'rel', 2, 2};
qnorms = {NaN, NaN, 'log'};
pqnames = {'Relative', 'Square', 'Regularized log-square'};
optimizers = {'somat3a', 'derandinfty', 'cmaes', 'ampso', 'apgskimode', 'madDE', 'fmincon'};
%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only NaN norms tried');
plen = length(pnorms);
qlen = length(qnorms);
olen = length(optimizers);
models = cell(N,olen,plen);
for i = 1:N
    dataN(i) = CreateBenchmarkProblem(rndErr, sysErr, minErr, id);
end

% For better parallelization we flatten three nested loops into a single
% one by (1) obtaining all the indices, (2) running computations in a
% flattened loop, (3) rewriting results to the original data structure
k = 1;
mk = cell(plen * olen *N, 1);
ijrep = [];
for i = 1:plen
    for j = 1:olen
        for rep = 1:N
            ijrep = [ijrep; i j rep k];
            k = k+1;
        end
    end
end

parfor k = 1:(plen * olen *N)
    i = ijrep(k, 1);
    j = ijrep(k, 2);
    rep = ijrep(k, 3);
    
    lambda_reg = [1 0];
    if(strcmp(qnorms{i},'log'))
        lambda_reg = [1, lambdas(4,2) * lossFunctionOrderMultiplier(4,2)];
    end
    
    optimizer = optimizers{j};
    
    data = dataN(rep);
    p = pnorms{i};
    q = qnorms{i};
    options = struct();
    
    try
        mk{k} = EstimateKineticModel(data,p,q,'batch',...
            lambda_reg, options, optimizer);
    catch
        warning('runBenchmarkingExperiment:Problem using function. Assigning a value of NaN.');
        mk{k} = NaN;
    end
    
end

k = 1;
ijrep = [];
for i = 1:plen
    for j = 1:olen
        for rep = 1:N
            models{rep,i,j} = mk{k};
            k = k+1;
        end
    end
end

save(savefilename);
disp(['Saved as ' savefilename])
end


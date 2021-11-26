function [xMin, fMin, flag, output] = DeRandInfty(fun, initialPopulation,...
    lBounds, uBounds, opts, varargin)
% DERANDINFTY - implementation of a variant of Differential Evolution
% algorithm called DE/rand/infty for real-parameter minimization with
% box constraints.
%
% Call:
% [xmin, fmin, flag, output] = DeRandInfty(fun, initialPopulation,...
%    lBounds, uBounds, opts)
%
% Input values:
% FUN - function to minimize
% INITIALPOPULATION - matrix containing initial population; each individual
% is encoded as a column vector
% LBOUNDS, UBOUNDS - vectors containing lower and upper bounds of the
% search space
% OPTS - structure with additional algorithm options  opts = struct('Dim',
% [problem dimension], 'MaxFunEvals', [maximal number of fitness function
% evaluation to use as stopping criterion], 'Disp', ['on'|'off' - displaying
% the intermediate results of computation])
%
% Output values:
% XMIN - argument, for which minimal value of function was found
% FMIN - minimal value of function found
% FLAG - currently not used (reserved for future functionality)
% OUTPUT - additional output information (e.g. FEs is the number of fitness
% function evaluations)
%
% Example 1:
% fun = @(x) sum(x.^2); % sphere function (x is a vector)
% dim = 10; % 10-dimensional optimization
% lBounds = -5*ones(dim,1); uBounds = 5*ones(dim,1); % optimization in hypercube [-5; 5]^10
% opts = struct('Dim', dim, 'MaxFunEvals', 1e4*dim);
% [xmin, fmin, flag, output] = DeRandInfty(fun, [], lBounds, uBounds, opts)
%
% Example 2:
% fun = @(x) sum(x.^2); % sphere function (x is a vector)
% dim = 30; % 30-dimensional optimization
% lBounds = -5*ones(dim,1); uBounds = 5*ones(dim,1); % optimization in hypercube [-5; 5]^30
% NP = 100*dim;
% initialPopulation = 2*rand(dim, NP)-1; % initialize population in hypercube [-1; 1]^30
% opts = struct('Dim', dim, 'MaxFunEvals', 1e4*dim);
% [xmin, fmin, flag, output] = DeRandInfty(fun, initialPopulation, lBounds, uBounds, opts)
%
% Description and discussion of the DE/rand/infty algorithm can be found
% in the following papers:
% 1. Opara K., Arabas J., Differential mutation based on population covariance
%    matrix, Lecture Notes in Computer Science 6238, pp. 114-123, Heidelberg, 2010
% 2. Opara K., Arabas J., Decomposition and Metaoptimization of Mutation
%    Operator in Differential Evolution, Lecture Notes in Computer Science,
%    Volume 7269/2012, pp. 110-118, DOI: 10.1007/978-3-642-29353-5_13, 2012
%
% Karol Opara, 26th July 2012

if (isempty(initialPopulation))
    dim = opts.Dim;
    if (isfield(opts,'PopSize') == false)
        popSize = 5*opts.Dim;
    else
        popSize = opts.PopSize;
    end
    pop = InitializePopulation.Uniform(lBounds,uBounds,popSize);
else
    pop = initialPopulation;
    [dim, ~] = size(pop);
end
if (isfield(opts,'MaxFunEvals'))
    maxFunEvals = opts.MaxFunEvals;
else
    maxFunEvals = 1e4*dim;
end
disp = true;
if (isfield(opts,'Disp') && strcmp(opts.Disp,'off'))
    disp = false;
end
fun = fcnchk(fun,length(varargin));

F = 0.7;
Cr = 0.9;

[fVal, funEvals] = evaluateObjectiveFunction(fun, pop, 0, varargin{:});

iter = 0;
while (true)
    [stop, flag] = isStopConditionMet(funEvals, maxFunEvals);
    if (stop == true)
        break;
    end
    
    % Variation operators
    % mutatnts = Mutation.DERand1(pop, F);
    mutatnts = Mutation.DeRandInfty(pop, F);
    crossedOverMutants = Crossover.Binomial(pop, mutatnts, Cr);
    
    % Constraint handling
    offsprings = Constraints.Reflection(crossedOverMutants, lBounds, uBounds);
   
    % Objective function evaluation
    [fValOffspring, funEvals] = evaluateObjectiveFunction(fun, offsprings,...
        funEvals, varargin{:});

    % Selection between parent and offspring
    offspringBetter = fValOffspring <= fVal;
    pop(:,offspringBetter) = offsprings(:,offspringBetter);
    fVal(offspringBetter) = fValOffspring(offspringBetter);
    
    iter = nextIteration(iter,disp,pop,fVal,funEvals);
end

[xMin, fMin, output] = prepareOutput(pop, fVal, funEvals);
end

function iter = nextIteration(iter,disp,pop,fVal,funEvals)
iter = iter + 1;
if (disp == false)
    return;
end
% reportIter = round(logspace(0,4,10));
reportIter = [1 2 5 1e1 3e1 5e1 1e2 3e2 5e2 1e3 3e3 5e3 1e4 3e4 5e4 1e5 ...
    3e5 5e5 1e6 3e6 5e6 1e7 3e7 5e7 1e8 3e8 5e8 1e9 3e9 5e9 1e10 3e10 5e10];
if (any(iter == reportIter))
    it = num2str(iter);
    fes = num2str(funEvals);
    fBest = num2str(min(fVal));
    cp = cov(pop.');
    popNorm = num2str(norm(cp));
    [~, d] = eig(cp);
    d = sort(diag(d),'descend');
    axRatio = abs(d(1)/d(2));
    axRat = num2str(axRatio);
    
    fprintf(['Iter = ' it ',\tFunEvals = ' fes ...
        ',\tFVal = ' fBest ',\tPopSpread = ' popNorm ',\t AxRatio = ' axRat '.\n']);
end
end

function [xMin, fMin, output] = prepareOutput(pop, fVal, funEvals)
[fMin, indMin] = min(fVal);
xMin = pop(:,indMin);
output.FunEvals = funEvals;
end

function [fVal, funEvals] = evaluateObjectiveFunction(fun, pop, funEvals, varargin)
[~, popSize] = size(pop);
funEvals = funEvals+popSize;
fVal = Inf(1,popSize);
parfor i = 1:popSize
    fVal(i) = fun(pop(:,i),varargin{:});
end
end

function [stop, flag] = isStopConditionMet(funEvals, maxFunEvals)
stop = false;
if (funEvals >= maxFunEvals)
    stop = true;
end
flag = ['Maximum number of function evaluations (' num2str(maxFunEvals) ') reached'];
end



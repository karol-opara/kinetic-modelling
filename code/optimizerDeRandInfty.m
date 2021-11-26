function [xMin, fMin, flag, output] = optimizerDeRandInfty(fun, initialPopulation,...
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
% License GNU GPL2+
%
% Karol Opara, 26th July 2012

if (isempty(initialPopulation))
    dim = opts.Dim;
    if (isfield(opts,'PopSize') == false)
        popSize = 5*opts.Dim;
    else
        popSize = opts.PopSize;
    end
    pop = InitializeUniform(lBounds,uBounds,popSize);
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
    mutatnts = MutationDeRandInfty(pop, F);
    crossedOverMutants = CrossoverBinomial(pop, mutatnts, Cr);
    
    % Constraint handling
    offsprings = ConstraintsReflection(crossedOverMutants, lBounds, uBounds);
    
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
for i = 1:popSize
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


function initialPopulation = InitializeUniform(lBounds, uBounds, popSize, randState)
% Uniform random initialization within given limits.
% Returns population matrix with each column denoting a
% single individual.

if nargin == 4
    RandStream.getDefaultStream.State = randState;
elseif nargin ~= 3
    error('Wrong number of arguments');
end

MinimizationProblemCheckBounds(lBounds, uBounds);
[dim, ~] = size(lBounds);

% calculate and return initial population matrix
initialPopulation = (repmat(lBounds, 1, popSize) + ...
    repmat((uBounds - lBounds), 1, popSize).*rand(dim, popSize));

end % initPopulation

function pop = MutationDeRandInfty(pop, F)
[~,popSize] = size(pop);
donors = pop(:,randi(popSize,1,popSize));
covMtx = 2*F^2*cov(pop.');
pop = MutationGaussian(donors,covMtx);
end

function pop = MutationGaussian(pop, covMtx)
% GAUSSIANMUTATION mutates each individual (column) of population pop
% by adding to it a realization of a normally distributed variable with
% zero mean and covariance matrix covMtx.
% Alternatively sigma^2 (square of std) can be given as input.

[dim, popSize] = size(pop);
[r, c] = size(covMtx);
% If covMtx is scalar make it a vector
if (r == 1 && c == 1)
    covMtx = covMtx * ones(dim,1);
end
mu = zeros(dim,1);
pop = pop + mvnrnd(mu, covMtx, popSize).';
end % Gaussian

function pop = CrossoverBinomial(pop, pop2, crossoverProbability)
% CROSSBINOMIAL Crossover two populations, substitute genes of individuals
% from base population with individuals
[dim, popSize] = size(pop);
cross = rand(dim,popSize) < crossoverProbability;
pop(cross) = pop2(cross);
end % crossBinomial

function pop = ConstraintsReflection(pop, lBounds, uBounds)
% REFLECTION hanles constraint breach by reflection.
% If single reflection is not enough, variable violating the
% constraint is set randomly within feasible set.

if(numel(lBounds) == 1 && numel(uBounds) == 1 ...
        && lBounds == -Inf && uBounds == Inf)
    return;
end

% If bounds are infinite, return.
if any(~isfinite(lBounds) | ~isfinite(uBounds))
    if any(isfinite(lBounds) | isfinite(uBounds))
        warning('ConstraintsReflection:constrReflection',...
            'Mix of finite and infinite bounds is not supported, no bounds assumed');
    end
    return;
end

[dim, popSize] = size(pop);
L = repmat(lBounds, 1, popSize);
U = repmat(uBounds, 1, popSize);
%             lowerConstr = pop < L;
%             pop = pop + lowerConstr.*(2*L - 2*pop);
%             upperConstr = pop > U;
%             pop = pop + upperConstr .* (2*U - 2*pop);
lowerConstr = pop < L;
pop(lowerConstr) = pop(lowerConstr) + (2*L(lowerConstr) - 2*pop(lowerConstr));
upperConstr = pop > U;
pop(upperConstr) = pop(upperConstr) + (2*U(upperConstr) - 2*pop(upperConstr));

% checking in case of large violations of constraints
largeViolation = pop < L | pop > U;
if (any(any(largeViolation)))
    nonFeasibleRatio = sum(sum(largeViolation))/(dim*popSize);
    if(nonFeasibleRatio > 0.01)
        warning('ConstraintsReflection:LargeViolation', ...
            'Large violations of constraints, uniform random resampling used instead');
    end
    randomPopulation = InitializeUniform(lBounds, uBounds, popSize);
    pop(largeViolation)= randomPopulation(largeViolation);
end
end % Reflection

function MinimizationProblemCheckBounds(lBounds, uBounds)
[rl cl] = size(lBounds);
[ru cu] = size(uBounds);
if (rl ~= ru)
    error('Lower and upper bounds must have the same dimension');
end
if(cl ~= 1 || cu ~= 1)
    error('Upper and lower bounds must be given as column vectors');
end
if (any(any(lBounds > uBounds)))
    error('Lower bounds cannot be greater than upper ones');
end
end
        
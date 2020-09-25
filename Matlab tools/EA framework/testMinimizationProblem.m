function test_suite = testMinimizationProblem
initTestSuite;
end

function testSetDefaultOptions
fitFun = @(x) sum(x);
initLBounds = [2 3 4].';
initUBounds = [3 4 5].';
minProb = MinimizationProblem(fitFun, initLBounds, initUBounds);
% Method SetDefaultOptions is called in constructor, therefore we can start
% asserting
assertTrue(minProb.EvalParallel);
assertEqual(-Inf,minProb.LBounds);
assertEqual(Inf,minProb.UBounds);
end

function testConstructor
FitnessFun = @(x) mean(x.^2)+1;
InitLBounds = [0 0 0].';
InitUBounds = [1 1 1].';

pOpts.Name = 'Sphere';
pOpts.LBounds = -Inf;
pOpts.UBounds = Inf;

minProb = MinimizationProblem(FitnessFun, InitLBounds, InitUBounds, pOpts);

assertEqual(minProb.Dimension, 3);
assertEqual(minProb.FitnessFunction([0 1; 0 1; 0 1]), [1 2]);
assertEqual(minProb.Name, 'Sphere');
assertEqual(minProb.LBounds, -Inf);
assertEqual(minProb.UBounds, Inf);
end

function testEvaluateFitnessFunction_Varargin_NoEvalParallel
fitFun = @(x,y) mean(x) + mean(y);
initLBounds = [2 3 4].';
initUBounds = [3 4 5].';
minProb = MinimizationProblem(fitFun, initLBounds, initUBounds);
minProb.EvalParallel = false;

pop = [1 1 1; 2 2 2; 3 3 3].';
% k-th column contains parameters for k-th individual from population pop
y = {10*pop(:,1), 10*pop(:,2), 10*pop(:,3)}; 

[fVal] = minProb.EvaluateFitnessFunction(pop, y);
assertEqual(fVal, [11 22 33]);
end

function testEvaluateFitnessFunction_Varargin_EvalParallel
fitFun = @(x,y) mean(x) + mean(y);
initLBounds = [2 3 4].';
initUBounds = [3 4 5].';
minProb = MinimizationProblem(fitFun, initLBounds, initUBounds);
minProb.EvalParallel = true;
pop = [1 1 1; 2 2 2; 3 3 3].';
y = 10*pop;
[fVal] = minProb.EvaluateFitnessFunction(pop, y);
assertEqual(fVal, [11 22 33]);
end

function testEvaluateFitnessFunction_EvalNonParallel
fitFun = @(x) mean(x.^2);
pop = [ 2 3 4 5 ];
pOpts.EvalParallel = false;

initLBounds = [2 3 4].';
initUBounds = [3 4 5].';
minProb = MinimizationProblem(fitFun, initLBounds, initUBounds, pOpts);

[fVal fMin xMin funEvals] = minProb.EvaluateFitnessFunction(pop);

assertEqual(fVal, [4 9 16 25]);
assertEqual(fMin, 4);
assertEqual(xMin, [2]);
assertEqual(funEvals, 4);
end

function testEvaluateFitnessFunction_EvalParallel
fitFun = @(x) sum(x.^2);
pop = [ 2 3 4 5
        2 3 4 5];
pOpts.EvalParallel = true;

initLBounds = [2 3 4].';
initUBounds = [3 4 5].';
minProb = MinimizationProblem(fitFun, initLBounds, initUBounds, pOpts);

[fVal fMin xMin funEvals] = minProb.EvaluateFitnessFunction(pop);

assertEqual(fVal, 2*[4 9 16 25]);
assertEqual(funEvals, 4);
assertEqual(xMin, [2 2].');
assertEqual(fMin, 8);
end
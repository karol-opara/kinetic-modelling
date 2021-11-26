function test_suite = testMinimizer
initTestSuite;
end
% 
% function testCheckArgs_OK
% checkArgs = getHandleByName('checkArgs');
% fitFun = @(x) x.^2;
% initRange = [1 2
%             3 5
%             3 9];
% [ok, dim, errMsg] = checkArgs(fitFun, initRange, []);
% assertEqual(dim, 3);
% assertTrue(ok);
% end
% 
% function testCheckArgs_WrongMatrixSize
% checkArgs = getHandleByName('checkArgs');
% fitFun = @(x) x.^2;
% initRange = [ 14 5 6
%             3 5 9];    
% [ok, dim, errMsg] = checkArgs(fitFun, initRange, []);
% assertFalse(ok)    
% end
% 
% function testCheckArgs_WrongInitRange
% checkArgs = getHandleByName('checkArgs');
% fitFun = @(x) x.^2;
% initRange = [ 14 6 % x1min =  14 < 6 x1max, hence error 
%             3 4];    
% [ok, dim, errMsg] = checkArgs(fitFun, initRange, []);
% assertFalse(ok)   
% end
% 
% function testCheckArgs_InfInitRange
% checkArgs = getHandleByName('checkArgs');
% fitFun = @(x) x.^2;
% initRange = [ -Inf 6
%             3 4];    
% [ok, dim, errMsg] = checkArgs(fitFun, initRange, []);
% assertFalse(ok)  
% end

% function h = getHandleByName(functionName)
% handles = Minimizer.GetSubfunctionHandles;
% h = NaN;
% fIndex = 0;
% for i=1:length(handles)
%     fName = func2str(handles{i});
%     if (strcmp(functionName, fName))
%         fIndex = i;
%         break;
%     end
% end
% assertFalse(fIndex == 0, 'The given function name not found in structure retuned by DMEA');
% h = handles{fIndex};
% end

function testSetUserOptions_MisspellingStruct
m = testMinimizerMock;
minProb = MinimizationProblem(@(x) sum(x), [2 3].', [4 5].');
m.SetDefaultOptions(minProb);

maxIter = m.StopCond.MaxIter;
warning('off', 'Minimizer:SetUserOptions:Structure');
inOpts.MISSPELLING_STRUCT.MaxIter = 44;
assertTrue(m.SetUserOptions(inOpts));
assertEqual(maxIter, m.StopCond.MaxIter);
warning('on', 'Minimizer:SetUserOptions:Structure');
end

function testSetUserOptions_OK
m = testMinimizerMock;
minProb = MinimizationProblem(@(x) sum(x), [2 3].', [4 5].');
m.SetDefaultOptions(minProb);
inOpts.StopCond.MaxIter = 44;
inOpts.AlgOpts.Name = '''unit test name''';
inOpts.Output.StopReason = '''ctr+c''';
warnings = m.SetUserOptions(inOpts);
assertEqual(m.StopCond.MaxIter, 44);
assertEqual(m.AlgOpts.Name, '''unit test name''');
assertEqual(m.Output.StopReason, '''ctr+c''');
assertEqual(warnings, false);
end

function testSetDefaultOptions
m = testMinimizerMock;
minProb = MinimizationProblem(@(x) sum(x), [1 2 3].', [4 5 6].');
m.SetDefaultOptions(minProb);
assertEqual(m.StopCond.MaxFunEvals,30000);
assertEqual(m.StopCond.MaxIter, 30000);
assertEqual(m.Output.Iter, 0);
assertEqual(m.Output.FMin, Inf);
assertEqual(m.AlgOpts.Name, 'Test algorithm');
assertEqual(m.AlgOpts.Version, '0.0.0.0');
end

function testUpdateOutput_LowerXMin
m = testMinimizerMock;

output.FMin = 3;
output.XMin = [4.1; 6.5];
output.FunEvals = 13;
output.Iter = 45;

m.Output = output;

funEvals = 4;
xMin = [2; 5];
fMin = 2;

m.UpdateOutput(funEvals, fMin, xMin);

output = m.Output;
assertEqual(output.FMin, 2);
assertEqual(output.XMin, [2; 5]);
assertEqual(output.FunEvals, 17);
assertEqual(output.Iter, 46);
end


function testUpdateOutput_GreaterXMin
m = testMinimizerMock;

m.Output.FMin = 1.6;
m.Output.XMin = [4.1; 6.5];
m.Output.FunEvals = 13;
m.Output.Iter = 0;
funEvals = 1;
xMin = [2; 5];
fMin = 2;

m.UpdateOutput(funEvals, fMin, xMin);

assertEqual(m.Output.FMin, 1.6);
assertEqual(m.Output.XMin, [4.1; 6.5]);
assertEqual(m.Output.FunEvals, 14);
assertEqual(m.Output.Iter, 1);
end


function testIsStopConditionMet_MaxFunEvals
m = testMinimizerMock;
m.StopCond.MaxFunEvals = 135;
m.Output.FunEvals = 137;
m.Output.Iter = 14;
m.StopCond.MaxIter = 100;
assertTrue(m.IsStopConditionMet());
assertEqual(m.Output.StopReason, 'MaxFunEvals exceeded');

m.StopCond.MaxFunEvals = 1345;
m.Output.FunEvals = 1345;
assertTrue(m.IsStopConditionMet());

m.StopCond.MaxFunEvals = 1345;
m.Output.FunEvals = 14;
assertFalse(m.IsStopConditionMet());
end

function testIsStopConditionMet_MaxIter
m = testMinimizerMock;
m.StopCond.MaxIter = 135;
m.Output.Iter = 137;
m.Output.FunEvals = 3;
m.StopCond.MaxFunEvals = 100;
assertTrue(m.IsStopConditionMet());
assertEqual(m.Output.StopReason, 'MaxIter exceeded');

m.StopCond.MaxIter = 1345;
m.Output.Iter = 1345;
assertTrue(m.IsStopConditionMet());

m.StopCond.MaxIter = 1345;
m.Output.Iter = 144;
assertFalse(m.IsStopConditionMet());
end

function testSetInitialPopulation
m = testMinimizerMock;
m.InitializeMinimizer();
initPop = [1.9 3.3; 3 2; -0.3 2.77].';
m.SetInitialPopulation(initPop);
pop = m.Pop;
assertEqual(initPop, pop(:,1:3));
end

% function testSpellingCheckInInOpts
% printDefaultOptions = getHandleByName('printDefaultOptions');
% getOptions = getHandleByName('getOptions');
% defOpts = printDefaultOptions();
% N = 10;
% inOpts.SpellingMistake = 123;
% dOpts = getOptions(inOpts, defOpts, N);
%
% dOptFields = fieldnames(dOpts);
% for i=1:length(dOptFields)
%     fieldName = dOptFields{i};
%     assertFalse(strcmp('SpellingMistake',fieldName));
% end
% end
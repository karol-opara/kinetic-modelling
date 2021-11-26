function test_suite = testUtils
initTestSuite;
end

function testSetUserOptions_MisspellingField
warning('off', 'Utils:SetUserOptions:Field');
userOpts.MISSSPELLING_FIELD = '''something''';
defOpts.Name = 'Def. opts name';
[opts, warnings] = Utils.SetUserOptions(defOpts, userOpts);
assertTrue(warnings);
assertEqual(opts, defOpts);
warning('on', 'Utils:SetUserOptions:Field');
end

function testSetUserOptions_OK
defOpts.A = 'Aa';
defOpts.Eleven = 11;
defOpts.NotANum = NaN;
defOpts.NotChanged = 342;
inOpts.Eleven = 44;
inOpts.A = 'Bb';
inOpts.NotANum = '''ctr+c''';
[opts warnings] = Utils.SetUserOptions(defOpts,inOpts);
assertFalse(warnings);
assertEqual(opts.Eleven, 44);
assertEqual(opts.A, 'Bb');
assertEqual(opts.NotANum, '''ctr+c''');
assertEqual(opts.NotChanged, 342);
end

function testParseDefaultOptions
optStr.Name =     '''DMEA''';
optStr.Version =  '''1.0.2.0''';
optStr.LastModificationDate = '''6-Sep-2011''';
optStr.PopSize =           '5*dim  % aka Np or mu';
optStr.ScalingFactor =     '0.6    % aka F';
optStr.EvaluatePopMidpoint='true   % evaluates a point at the middle of population';

minProb.Dimension = 3;

opts = Utils.ParseDefaultOptions(optStr, minProb);

dOpts.Name = 'DMEA';
dOpts.Version = '1.0.2.0';
dOpts.LastModificationDate = '6-Sep-2011';
dOpts.PopSize = 15;
dOpts.ScalingFactor = 0.6;
dOpts.EvaluatePopMidpoint = true;

assertEqual(opts, dOpts);
end
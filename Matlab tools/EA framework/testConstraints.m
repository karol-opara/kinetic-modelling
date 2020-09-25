function test_suite = testConstraints
initTestSuite;
end

function testReflection_NoConstraints
pop = rand(3,4);
lBounds = -Inf;
uBounds = Inf;
pop1 = Constraints.Reflection(pop, lBounds, uBounds);
assertEqual(pop, pop1);
end

function testReflection_LBound
pop =       [ 1 3 5 7
              7 1 9 6];
lBounds1 = ones(2,1);
uBounds1 = ones(2,1)*9;

pop1 = Constraints.Reflection(pop, lBounds1, uBounds1);
assertEqual(pop, pop1); 
end

function testReflection_BothBounds
pop =       [ 1 3 5 7
              7 1 9 6];
lBounds2 = ones(2,1)*2;
uBounds2 = ones(2,1)*6;
popExpected = [3 3 5 5
               5 3 3 6];
           
pop2 = Constraints.Reflection(pop, lBounds2, uBounds2);
assertEqual(popExpected, pop2); 
end

function testReflection_LargeViolation
popLargeViolation = [2.2  3e4+0.4  0.4
                     1   33  0.4];
lBounds3 = [0 0].';
uBounds3 = [1 1].';

warning('off', 'ConstraintsReflection:LargeViolation');
pop3 = Constraints.Reflection(popLargeViolation, lBounds3, uBounds3);
warning('on', 'ConstraintsReflection:LargeViolation');
popExpected3 = [0.2 0.4 0.4
                1 1 0.4];
            
% assert all individuals are within feasible set
assertTrue(all(all(pop3 <= ones(2,3) & pop3 >= zeros(2,3))), 'Not all individuals are within feasible set');
            
% assert violations are handled only when necessary
assertEqual([pop3(1,3) pop3(2,1) pop3(2,3)],...
    [popLargeViolation(1,3) popLargeViolation(2,1) popLargeViolation(2,3)]);
end

function testReflection_1D
pop4 = [1 2 3 4 5 6 7];
lBounds4 = [2.5];
uBounds4 = [7];
pop4res = Constraints.Reflection(pop4, lBounds4, uBounds4);
pop4expected = [4 3 3 4 5 6 7];
assertEqual(pop4res, pop4expected);
end
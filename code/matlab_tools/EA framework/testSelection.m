function test_suite = testSelection
initTestSuite;
end

function testProportional_proportionality
popSize = 3*4e3;
parentPop = repmat([1 2 3; 4 5 6; 7 8 9].',1,popSize/3);
parentFitness = repmat([-1 0 1], 1, popSize/3);
offspring = Selection.Proportional(parentPop, parentFitness);

offspring123Count = sum(sum(repmat([1 2 3].',1,popSize) == offspring));
assertEqual(offspring123Count,0);

offspring456Count = sum(sum(repmat([4 5 6].',1,popSize) == offspring));
assertAlmostEqual(offspring456Count/3, popSize/3, 0.15);

offspring789Count = sum(sum(repmat([7 8 9].',1,popSize) == offspring));
assertAlmostEqual(offspring789Count/3, 2* popSize/3, 0.15);
end

function testProportional_offspringNumber
popSize = 234;
parentPop = repmat([1 2 3; 4 5 6; 7 8 9].',1,popSize/3);
parentFitness = repmat([-1 0 1], 1, popSize/3);
offspring = Selection.Proportional(parentPop, parentFitness, 345);
[~, offSize] = size(offspring);
assertEqual(offSize, 345);
end

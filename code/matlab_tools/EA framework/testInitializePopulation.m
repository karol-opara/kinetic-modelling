function test_suite = testInitializePopulation
initTestSuite;
end

function testUniform_Multidimensional
initRange = [-1      5
              0      7
              1     11
              2     14
              -12.4  12];
lBounds = initRange(:,1);
uBounds = initRange(:,2);
popSize = 1e4;       
  
pop = InitializePopulation.Uniform(lBounds, uBounds, popSize);

lowerBoundMatrix = repmat(initRange(:,1),1,popSize);
upperBoundMatrix = repmat(initRange(:,2),1,popSize);

assertFalse(any(any(pop < lowerBoundMatrix)));
assertFalse(any(any(pop > upperBoundMatrix)));

for i = 1:10
    decileCount(i) = sum(pop(3,:)>=i & pop(3,:)<i+1);
end
assertAlmostEqual(decileCount, ones(1,10)*popSize/10, 0.15);
end

function testUniform_Onedimensional
initRange = [-1      9];
popSize = 1e4;       
  
pop = InitializePopulation.Uniform(-1, 9, popSize);

lowerBoundMatrix = repmat(-1,1,popSize);
upperBoundMatrix = repmat(9,1,popSize);

assertFalse(any(any(pop < lowerBoundMatrix)));
assertFalse(any(any(pop > upperBoundMatrix)));

for i = -1:8
    decileCount(i+2) = sum(pop>=i & pop<i+1);
end
assertAlmostEqual(decileCount,  ones(1,10)*popSize/10, 0.15);
end



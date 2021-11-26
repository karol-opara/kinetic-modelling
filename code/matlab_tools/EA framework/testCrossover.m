function test_suite = testCrossover
initTestSuite;
end

function testBinomial_FixedRandStream
popBase = [ 1 3 5 7
            7 1 9 3];
popTarget = [ 2 4 8 0
              6 4 6 2];
crossoverProbability = 0.3;
            
defaultStream = RandStream.getDefaultStream;
          
stream = RandStream('mrg32k3a');
RandStream.setDefaultStream(stream);
% within this substream we have
% >> rand(2,3)
% >> ans =
%     0.7270    0.9387    0.0277    0.2112
%     0.4522    0.2360    0.1316    0.9324

pop = Crossover.Binomial(popBase, popTarget, crossoverProbability);
% get back to defauld stream
RandStream.setDefaultStream(defaultStream);
     
predictedOffspring = [ 1 3 8 0
                       7 4 6 3]; 
assertEqual(predictedOffspring, pop);

end
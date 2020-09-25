function test_suite = testMutation
initTestSuite;
end

function testDeRandInfty_1D
       mu = [1 -1]; Sigma = [.9 .4; .4 .3];
       pop = mvnrnd(mu, Sigma, 500).';
       F = 3.4;
       
       pop2 = Mutation.DeRandInfty(pop, F);
       %plot(pop(1,:),pop(2,:),'ro',pop2(1,:),pop2(2,:),'b.');
       
       assertElementsAlmostEqual(mean(pop2.'),mu,'relative',0.1);
       assertElementsAlmostEqual(cov(pop2.'),2*F*F*Sigma, 'relative', 0.1);
end

function testDERand1_1D
pop1D = [1 3 7];
F = 10;
pop1DRes = Mutation.DERand1(pop1D, F);
pop1DResA = [-39 -57 -13];
pop1DResB = [41 63 27];
assertTrue(all(all(pop1DRes == pop1DResA | pop1DRes == pop1DResB)));
end

function testDERand1_2D
pop2D = [1 3 7
         1 3 7];
F = 10;
pop2DRes = Mutation.DERand1(pop2D, F);
pop2DResA = [-39 -57 -13
             -39 -57 -13];
pop2DResB = [41 63 27
             41 63 27];
assertTrue(all(all(pop2DRes == pop2DResA | pop2DRes == pop2DResB)));
end

function testDERand1_largePopulation
popLarge = rand(10,1e3);
popLargeRes = Mutation.DERand1(popLarge,0.9);
assertTrue(all(all(popLarge ~= popLargeRes)));
end

function testGaussian_1D
pop = zeros(1,1e4);
sigma = 33;
popOffspring = Mutation.Gaussian(pop, sigma);
pdn = ProbDistUnivParam('normal',[0 sqrt(sigma)]);
[h p]= kstest(popOffspring,pdn,1e-3,'unequal'); %Kolmogorov-Smirnoff test
assertEqual(h,false); % assert distributions are not unequal
end

function testGaussian_2D
pop = zeros(2,1e4);
sigma = [9 4; 4 3];
popOffspring = Mutation.Gaussian(pop, sigma);
[V D] =eig(cov(popOffspring.'));
assertAlmostEqual(V,[0.45 -0.9; -0.9 -0.45],0.15);
assertAlmostEqual(diag(D),[1; 11], 0.15);
end
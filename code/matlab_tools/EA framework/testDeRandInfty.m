function test_suite = testDeRandInfty
initTestSuite;
end


function testDeRandInfty_SphereConvergence
fun = @(x) sum(x.^2); % sphere function (x is a vector)
dim = 5; % 10-dimensional optimization
lBounds = -5*ones(dim,1); uBounds = 5*ones(dim,1); % optimization in hypercube [-5; 5]^10
opts = struct('Dim', dim, 'MaxFunEvals', 1e3*dim, 'PopSize', 50,'Disp','off');
[xmin, fmin, flag, output] = DeRandInfty('tmp2', [], lBounds, uBounds, opts);
assert(fmin<0.01);
assert(norm(xmin)<10);
end


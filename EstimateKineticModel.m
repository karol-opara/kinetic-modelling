function model = EstimateKineticModel(data,p,q,type,lambda,options,weights)
if (nargin == 0)
    warning('EstimateKineticModel:NoInputArguments',...
        'No input arguments given, using the default ones')
    %     tmp = getExperimentalData();
    %     data = tmp{5};
    load saveExperimentalData20120807;
    data=data{5};
    p = 2;
    q = 2;
    type = 'batch';
elseif (nargin == 2)
    q = 2;
    warning('EstimateKineticModel:NoOuterNorm',...
        'Using default value q=2 for the outer norm');
elseif (nargin == 3)
    lambda = [1 1e-2 1e-1 1 1e-1 1e-1];
    warning('EstimateKineticModel:NoLagrangeMultipliers',...
        'No lagrange multipliers, using default ones');
end
if (nargin<6)
    options = struct();
end
if (nargin<7)
    weights = [1 1 1 1 1 1];
    %warning('EstimateKineticModel:NonuniformWeights',['No weights provided, using the default: ' num2str(weights)]);
end
%plotExprimentalData(data);

z0 = data.z0ml([1 3 4 5]).'; % MeOH, TG, DG, MG
k0 = rand(6,1)*0.1;

lBounds = [1e-8*ones(1,6)].';%, 5,  0.2, 1e-6, 1e-6].';
uBounds = [4.0E+00	2.7E+01	5.5E+01	6.6E+01	2.4E+00	1.8E-01].'; % the maximum from other papers

tic
% [k, z0opt] = deRandInftyOptimization(length(k0), lBounds, uBounds,...
%     data.zml, data.timeZ, data.z0ml,p,q);
[k, z0opt, out] = cmaesOptimization(length(k0), lBounds, uBounds,...
    data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
% [k, z0opt, out] = fMinConOptimization(length(k0), lBounds, uBounds,...
%     data.zml, data.timeZ, data.z0ml,p,q,type,lambda);
time = toc();

model = struct('data', data, 'k', k, 'z0opt', z0opt, 'type', type, 'executionTime',toc, 'optimizerOutput', out);
end

function [k, z0opt, fMinConOut] = fMinConOptimization(~, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda)
options=optimset('fmincon');
options.Algorithm = 'active-set';
[k,~,~,fMinConOut] = fmincon('ObjectiveFunction',z0ml,[],[],[],[],lBounds,uBounds,[],options,zml, timeZ, z0ml, p, q, type,lambda);
z0opt= z0ml;
end

function [k, z0opt] = deRandInftyOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type)
opts.Dim = dim;
%warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to very low value -- use for debug only');
opts.MaxFunEvals = 1e4*dim;
opts.PopSize = 5*dim;
opts.MinPopNorm = 1e-3;
[k, ~]= DeRandInfty('ObjectiveFunction', [], lBounds, uBounds, opts, zml, timeZ, z0ml, p, q,type);
z0opt = z0ml;
end

function [k, z0opt, cmaesOut] = cmaesOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts = cmaes('defaults');
opts.MaxFunEvals = 1e4*dim;
% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;
if(isempty(fieldnames(options))==false) % i.e. ther is options structure provided by user
    opts.MaxFunEvals = options.MaxFunEvals;
end
opts.MaxIter = Inf;
opts.LBounds = lBounds;
opts.UBounds = uBounds;
opts.SaveVariables = 'off';
opts.DispFinal = 0;
opts.DispModulo = 0;
opts.SaveFilename = 0;
opts.LogFilenamePrefix = 0;
opts.LogTime = 1;
warning('off','cmaes:logging');
opts.Restarts = 5;
[k, ~, ~, ~, cmaesOut, ~]= cmaes('ObjectiveFunction', (lBounds+uBounds)/2, [], opts, zml, timeZ, z0ml, p, q, type,lambda,weights);
z0opt = z0ml;
end
function model = EstimateKineticModel(data,p,q,type,lambda,options, optimizer, weights)
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
    optimizer = 'cmaes';
end
if (nargin<8)
    weights = [1 1 1 1 1 1];
    %warning('EstimateKineticModel:NonuniformWeights',['No weights provided, using the default: ' num2str(weights)]);
end
%plotExprimentalData(data);

z0 = data.z0ml([1 3 4 5]).'; % MeOH, TG, DG, MG
k0 = rand(6,1)*0.1;

lBounds = [1e-8*ones(1,6)].';%, 5,  0.2, 1e-6, 1e-6].';
uBounds = [4.0E+00	2.7E+01	5.5E+01	6.6E+01	2.4E+00	1.8E-01].'; % the maximum from other papers

tic
switch optimizer
    case 'madDE'
        [k, z0opt, out] = madDeOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    case 'fmincon'
        [k, z0opt, out] = fMinConOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda);
    case 'cmaes'
        [k, z0opt, out] = cmaesOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    case 'derandinfty'
        [k, z0opt, out] = deRandInftyOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    case 'apgskimode'
        [k, z0opt, out] = apgskimodeOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    case 'somat3a'
        [k, z0opt, out] = somat3aOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    case 'ampso'
        [k, z0opt, out] = ampsoOptimization(length(k0), lBounds, uBounds,...
            data.zml, data.timeZ, data.z0ml,p,q,type,lambda,options,weights);
    otherwise
        warning('EstimateKineticModel:EstimateKineticModel',...
            ['Unsupported optimizer type: ' optimizer]);
end

time = toc();

model = struct('data', data, 'k', k, 'z0opt', z0opt, 'type', type,...
    'executionTime',toc, 'optimizerOutput', out);
end

function [k, z0opt, fMinConOut] = fMinConOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda)
options=optimset('fmincon');
options.Algorithm = 'active-set';
options.MaxFunctionEvaluations = 1e4*dim;

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% options.MaxFunctionEvaluations = 1e1*dim;

[k,~,~,fMinConOut] = fmincon('ObjectiveFunction',z0ml,[],[],[],[],lBounds,uBounds,[],options,zml, timeZ, z0ml, p, q, type,lambda);
z0opt= z0ml;
end

function [k, z0opt, out] = deRandInftyOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts.Dim = dim;
opts.Disp = 'off';
opts.MaxFunEvals = 1e4*dim;
opts.PopSize = 5*dim;
opts.MinPopNorm = 1e-3;

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;

[k, ~, ~, out]= optimizerDeRandInfty('ObjectiveFunction', [], lBounds, uBounds, opts, zml, timeZ, z0ml, p, q,type,lambda,weights);
z0opt = z0ml;
end

function [k, z0opt, out] = somat3aOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
SOMApara.PopSize    = 1500;
SOMApara.N_jump     = 100;
SOMApara.m     		= 5;
SOMApara.n     		= 4;
SOMApara.k     		= 5;
Info.f_star         = 1.000000000;
Info.the_func       = 'ObjectiveFunction';
Info.FEs_Max        = 1e4*dim;
Info.dimension      = 6;
Info.Search_Range   = [-8192, 8192];

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;

out=NaN;
[Best , array_digit , FEs , Mig] = optimizerSOMA_T3A(Info, SOMApara, 'ObjectiveFunction', lBounds, uBounds, zml, timeZ, z0ml, p, q,type,lambda,weights);
k = Best.Positon.';
z0opt = z0ml;
end

function [k, z0opt, out] = apgskimodeOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts.Dim = dim;
opts.MaxFunEvals = 1e4*dim;
opts.PopSize = 5*dim;
opts.MinPopNorm = 1e-3;

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;

[k,Best_fit,All_fit,nfes,out]= optimizerAPGSK_IMODE('ObjectiveFunction', dim, opts.MaxFunEvals, lBounds, uBounds, zml, timeZ, z0ml, p, q,type,lambda,weights);
z0opt = z0ml;
end

function [k, z0opt, out] = madDeOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts.Dim = dim;
opts.MaxFunEvals = 1e4*dim;
opts.PopSize = 5*dim;
opts.MinPopNorm = 1e-3;

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e2*dim;

out=struct();
[k, best_val]= optimizerMadDE('ObjectiveFunction', dim, opts.MaxFunEvals, lBounds, uBounds, zml, timeZ, z0ml, p, q,type,lambda,weights);
out.best_val = best_val;
z0opt = z0ml;
end

function [k, z0opt, out] = ampsoOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts.Dim = dim;
opts.MaxFunEvals = 1e4*dim;
opts.PopSize = 5*dim;
opts.MinPopNorm = 1e-3;

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;
[k, best_fit, out] = optimizerAMPSO('ObjectiveFunction', opts.Dim, lBounds, uBounds, opts.MaxFunEvals,  zml, timeZ, z0ml, p, q, type, lambda, weights);
z0opt = z0ml;
end

function [k, z0opt, cmaesOut] = cmaesOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,options,weights)
opts = optimizerCMAES('defaults');
opts.MaxFunEvals = 1e4*dim;
if(isempty(fieldnames(options))==false) % i.e. there is the 'options' structure provided by user
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

% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;

[k, ~, ~, ~, cmaesOut, ~]= optimizerCMAES('ObjectiveFunction', (lBounds+uBounds)/2, [], opts, zml, timeZ, z0ml, p, q, type,lambda,weights);
z0opt = z0ml;
end
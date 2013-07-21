function KinteticsModelLandscape(name, condition, p, q, lambda)
% tmp = getExperimentalData();
% data = tmp{5};
load saveExperimentalData20120807
data = data{condition};

type = 'batch';

dir = 'Results/';
savefilename = [dir 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    '_KineticsModelLandscape_condition_' num2str(condition) '_' name '_L' num2str(p)...
    '_Q' num2str(q)];
 tstart = tic;


z0 = data.z0ml.';

N = 0;
ulimf = 5;
llimf = 0;

ulimb = 10;
llimb = 0;

linspf = linspace(llimf, ulimf, N);
lgspf = [2e-2 5e-2 1e-1 2e-1 5e-1 1e0 2e0 5e0 1e1 2e1 5e1 1e2 2e2];
spf = unique(sort([linspf lgspf]));

linspb = linspace(llimb, ulimb, N);
lgspb = [2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1 1e0 2e0 5e0 1e1 2e1 5e1 1e2 2e2];
spb = unique(sort([linspb lgspb]));

fwd = spf;
bck = spb;

z0 = data.z0ml;
zml = data.zml;
zTime = data.timeZ;

[kForward, kBackward] = meshgrid(fwd,bck);
kForward=kForward.';
kBackward=kBackward.';
for i=1:length(spf)
    parfor j=1:length(spb)
        sumForward = kForward(i,j);
        sumBackward = kBackward(i,j);
        %tic
        [K{i,j}, ErrLp(i,j), Flag{i,j}, Output{i,j}] = CalculateErrorsK...
            (sumForward, sumBackward, z0, zml, zTime,p, q, type, lambda); 
        %toc
        KProjected{i,j} = projectK(K{i,j}, sumForward, sumBackward);
    end
    fprintf('.');
    save([savefilename '_partial']);
end


telapsed = toc(tstart);
[h, m, s]=sec2hms(telapsed);
disp(['Elapsed time is ' num2str(h) ':' num2str(m) ':' num2str(s)]);
save(savefilename);
disp(['Saved optimization results as ' savefilename '.mat']);
end

function [k, err, flag, output] = CalculateErrorsK(sumForward, sumBackward, z0, zexp, timeZ,p,q,type,lambda)
lBounds = 1e-4*ones(6,1);
uBounds = 2.1e2*ones(6,1);
dim = 6; 
[k, err, flag, output] = cmaesOptimization(dim, lBounds, uBounds,...
    zexp, timeZ, z0,p,q,type,lambda,sumForward, sumBackward);
% [k, err, flag, output] = deRandInftyOptimizationK(6, lBounds, uBounds,...
%     z0, zexp, timeZ, sumForward, sumBackward,p);
end

function [k, err, flag, cmaesOut] = cmaesOptimization(dim, lBounds, uBounds, zml, timeZ, z0ml, p, q, type,lambda,sumForward, sumBackward)
opts = cmaes('defaults');
opts.MaxFunEvals = 1e4*dim;
% warning('EstimateKineticModel:MaxFunEvals','Max fun evals set to low value -- use for debug only');
% opts.MaxFunEvals = 1e1*dim;
opts.MaxIter = Inf;
opts.LBounds = lBounds;
opts.UBounds = uBounds;
opts.SaveVariables = 'off';
opts.DispFinal = 0;
opts.DispModulo = 0;
opts.SaveFilename = 0;
      0*log1p(norm(k,p))+...
	  0.01*sum(abs(k-sumForward-sumBackward));
warning('off','cmaes:logging');
opts.Restarts = 4;
%opts.PopSize = 24;
[k, err, ~, flag, cmaesOut, ~]= cmaes('KineticsModelLandscapeConstrainedObjectiveFunctionK',...
    (lBounds+uBounds)/2, [], opts, zml, timeZ, z0ml, p, q, type,lambda,sumForward, sumBackward);
end

% function [k, err, flag, output] = deRandInftyOptimizationK(dim, lBounds, uBounds, z0, zml,...
%     timeZ, sumForward, sumBackward,p)
% opts.Dim = dim;
% opts.MaxFunEvals = 2e5; %10h
% opts.MinPopNorm = 1e-4;
% opts.lBounds = lBounds;
% opts.uBounds = uBounds;
% opts.PopSize = 15*dim;%10*dim;
% [k, err, flag, output]= DeRandInfty(@ConstrainedObjectiveFunctionK, ...
%     [], lBounds, uBounds, opts, z0, zml, timeZ, sumForward, sumBackward,p);
% end

function k = projectK(k, sumForward, sumBackward)
ifw = [1 3 5];
ibk = [2 4 6];
k(ifw) = sumForward .* k(ifw)./sum(k(ifw));
k(ibk) = sumBackward .* k(ibk)./sum(k(ibk));
end
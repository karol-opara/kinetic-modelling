function KinteticsModelLandscape(name, condition, p)
% tmp = getExperimentalData();
% data = tmp{5};
load saveExperimentalData20120807
data = data{condition};

dir = 'Results/';
savefilename = [dir 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    '_KineticsModelLandscape_condition_' num2str(condition) '_' name '_L' num2str(p)];
 tstart = tic;


z0 = data.z0ml.';

N = 0;
ulimf = 5;
llimf = 0;

ulimb = 10;
llimb = 0;

linspf = linspace(llimf, ulimf, N);
lgspf = [2e-2 5e-2 1e-1 2e-1 5e-1 1e0 2e0 5e0 1e1 2e1 5e1 1e2];
spf = unique(sort([linspf lgspf]));

linspb = linspace(llimb, ulimb, N);
lgspb = [2e-2 5e-2 1e-1 2e-1 5e-1 1e0 2e0 5e0 1e1 2e1 5e1 1e2 2e2];
spb = unique(sort([linspb lgspb]));

fwd = spf;
bck = spb;

z0 = data.z0ml;
zml = data.zml;
zTime = data.timeZ;

i=1;
[kForward, kBackward] = meshgrid(fwd,bck);
kForward=kForward.';
kBackward=kBackward.';
for i=1:length(spf)
    parfor j=1:length(spb)
        sumForward = kForward(i,j);
        sumBackward = kBackward(i,j);
        tic
        [K{i,j}, ErrLp(i,j), Flag{i,j}, Output{i,j}] = CalculateErrorsK...
            (sumForward, sumBackward, z0, zml, zTime,p); 
        toc
        KProjected{i,j} = projectK(K{i,j}, sumForward, sumBackward);
    end
    save([savefilename '_partial']);
end


telapsed = toc(tstart);
[h, m, s]=sec2hms(telapsed);
disp(['Elapsed time is ' num2str(h) ':' num2str(m) ':' num2str(s)]);
save(savefilename);
disp(['Saved optimization results as ' savefilename '.mat']);
end

function [k, err, flag, output] = CalculateErrorsK(sumForward, sumBackward, z0, zexp, timeZ,p)
lBounds = 1e-4*ones(6,1);
uBounds = 1e2*ones(6,1);

[k, err, flag, output] = deRandInftyOptimizationK(6, lBounds, uBounds,...
    z0, zexp, timeZ, sumForward, sumBackward,p);
end

function [k, err, flag, output] = deRandInftyOptimizationK(dim, lBounds, uBounds, z0, zml,...
    timeZ, sumForward, sumBackward,p)
opts.Dim = dim;
opts.MaxFunEvals = 2e5; %10h
opts.MinPopNorm = 1e-4;
opts.lBounds = lBounds;
opts.uBounds = uBounds;
opts.PopSize = 15*dim;%10*dim;
[k, err, flag, output]= DeRandInfty(@ConstrainedObjectiveFunctionK, ...
    [], lBounds, uBounds, opts, z0, zml, timeZ, sumForward, sumBackward,p);
end

function err = ConstrainedObjectiveFunctionK(k, z0, zml, timeZ, sumForward, sumBackward,p)
k = projectK(k, sumForward, sumBackward);

[t, z] = SolveKineticModel(z0,k);
% HERE: compute fit error of the kinetic model based on experimental z's

zKinetic = interp1(t,z,timeZ);
%p = 0.5;%0.25;%0.5;

ferr = (sum(abs(zKinetic - zml).^p)).^(1/p);
fitErr = norm(ferr,2);
err = 1*(fitErr) +... % fit error
      0*log1p(norm(k,p))+...
	  0.01*sum(abs(k-sumForward-sumBackward));
end


% function [k, z0, err] = CalculateErrorsKz(sumForward, sumBackward, zexp, timeZ)
% lBounds = [1e-4*ones(1,6), 5,  0.2, 1e-6, 1e-6].';
% uBounds = [1e2*ones(1,6),	15, 0.8, 1e-1, 1e-2].';
% 
% [k, z0, err] = deRandInftyOptimization(10, lBounds, uBounds,...
%     zexp, timeZ, sumForward, sumBackward);
% end
% 
% function [k, z0opt, err] = deRandInftyOptimizationKz(dim, lBounds, uBounds, zml,...
%     timeZ, sumForward, sumBackward)
% opts.Dim = dim;
% opts.MaxFunEvals = 2.2e5;
% opts.lBounds = lBounds;
% opts.uBounds = uBounds;
% opts.PopSize = dim+3;%10*dim;
% [kz, fMin]= DeRandInfty(@ConstrainedObjectiveFunctionKz, ...
%     [], lBounds, uBounds, opts, zml, timeZ, sumForward, sumBackward);
% k = kz(1:6);
% z0opt = zeros(6,1);
% z0opt([1 3 4 5]) = kz(7:10);
% err = fMin;
% end
% 
% function err = ConstrainedObjectiveFunctionKz(kz, zml, timeZ, sumForward, sumBackward)
% %fprintf('.');
% k= kz(1:6);
% 
% k = projectK(k, sumForward, sumBackward);
% 
% z0opt = zeros(6,1);
% z0opt([1 3 4 5]) = kz(7:10);
% 
% [t, z] = SolveKineticModel(z0opt,k);
% % HERE: compute fit error of the kinetic model based on experimental x's
% 
% zKinetic = interp1(t,z,timeZ);
% p = 2;%0.25;%0.5;
% ferr = (sum(abs(zKinetic - zml).^p)).^(1/p);
% fitErr = ferr*[1 0 1 1 1 1].';
% 
% err = 1*(fitErr); % fit error
% %      1e2*(abs(sumForward-sum(k([1 3 5])))+abs(sumBackward-sum(k([2 4 6])))).^2;
% end

function k = projectK(k, sumForward, sumBackward)
ifw = [1 3 5];
ibk = [2 4 6];
k(ifw) = sumForward .* k(ifw)./sum(k(ifw));
k(ibk) = sumBackward .* k(ibk)./sum(k(ibk));
end
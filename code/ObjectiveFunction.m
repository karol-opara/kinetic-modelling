function err = ObjectiveFunction(k, zml, timeZ, z0opt,p,q,type, lambda, weights)
if nargin < 9
    weights = ones(1,6);
end
if nargin < 6
    q = 2;
    warning('ObjectiveFunction:NoOuterNorm',...
        'Using default value q=2 for the outer norm');
    type = 'batch';
    warning('ObjectiveFunction:NoType',...
        'Using default value type=''batch'' for the outer norm');
end
if nargin < 8
    warning('ObjectiveFunction:NoLagrangeMultipliers',...
        'No Lagrange multipliers, using the default ones');
    lambda = [1 1e-2 1e-1 1 1e-1 1e-1];
end

if strcmp(type,'membrane')
    timeSpan = [0 200];
else
    timeSpan = [0 100];
end

[t, z] = SolveKineticModel(z0opt,k,type,timeSpan);
% HERE: compute fit error of the kinetic model based on experimental x's

zKinetic = interp1(t,z,timeZ);
%p = 2;%0.25;%0.5;

%[zml, zKinetic] = studentizeData(zml, zKinetic);

errorMatrix = zKinetic-zml;
[r, c] = size(errorMatrix);
errorMatrix = repmat(weights,r,1).*errorMatrix;

if (isnan(q))
    if(ischar(p) && strcmp(p,'rel'))
        relativeErrorMatrix = abs(zKinetic-zml)./max(zml,1e-8); % we are replacing zeros with a small constant 1e-8 to enable division
        fitErr = sum(reshape(relativeErrorMatrix,r*c,1));
    else
        fitErr = myNorm(reshape(errorMatrix,r*c,1),p);
    end
else
    ferr = myNorm(errorMatrix,p);
    fitErr = myNorm(ferr,q);
end

% we want to minimize: fit error and fit in constraints
err = lambda(1)*fitErr +... % fit error
    lambda(2)*norm(k,2);%+... % reaction rates should not be too large
%lambda(2)*log1p(var(sum(zKinetic.'))) + ... % the total number of moles should not much (if volume is constant)
%     lambda(3)*log1p(sum(sum(abs(z.').*(z.'<0)))) +... % concentrations must be nonnegative
%     lambda(4)*log1p(sum(abs(k)))+... % reaction rates should not be too large
%     lambda(5)*log1p(sum(1./(1e-4+k))); % we do not want trivial solution k1=k2=...=k6=0
%lambda(4)*log1p(sum((x0opt(:,2)-x0opt(:,1)).*(x0opt(:,2)>x0opt(:,1)))); % we have more TG than MG

%err = fitErr;
if (rand()<1e-4)
    pause(1e-3);
end
end

function y = myNorm(x,p)
if (ischar(p))
    if(strcmp(p,'log'))
        y = sum(log1p(abs(x)));
    else
        error('ObjectiveFunction:myNorm','Inappropriate norm specified');
    end
else
    y = (sum(abs(x).^p)).^(1/p);
end
end

function [yExpS, yModS] = studentizeData(yExp,yMod)
[rs, ~] = size(yExp);
[rm, ~] = size(yMod);

yStd = std(yExp);
yMean = mean(yExp);

yExpS = (yExp-repmat(yMean,rs,1))./repmat(yStd,rs,1);
yModS = (yMod-repmat(yMean,rm,1))./repmat(yStd,rm,1);
end




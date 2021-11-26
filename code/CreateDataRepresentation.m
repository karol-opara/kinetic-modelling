
function dataRepresentation = CreateDataRepresentation(name)
% load Results\save_KineticModels_2012-08-16_172142_DataFrom20120807_L2_Lp
load('Results/save_KineticModels_2013-03-19_113532_8of13experiments_NRTLvalidation.mat')


tstart = tic;

t = 2:0.5:70;
N = 25;

savefilename = ['Results/' 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') name];

% k = getHandmadeKineticModelFit();
% warning('Handmade k values used');

%models = models(:,2);
for i=1:length(models)
    %d = data{i};
    d = models{i}.data;
    [tk, zk] = SolveKineticModel(models{i}.z0opt,models{i}.k);
    z = interp1(tk,zk,d.timeZ);
    zf = interp1(tk,zk,t);
    
    [alpha{i} alphaExp{i}] = ApproximateVolumes(d.xVolume, d.timeX, d.yVolume, d.timeY,t);
    %     [xf, yf] = FitCurves(data.xml, data.timeX, data.yml, data.timeY, z, data.timeZ, xVol, yVol, t);
    
    [xf, yf, yzmf] = SplitZToXY(d.xml,d.timeX, d.yml,d.timeY, zf, t, alpha{i},t);
    
    [xrepr, yrepr, tr] = TakeNewData(xf,yf,t,N);
    xrecon = interp1(t,xf,d.timeX);
    yrecon = interp1(t,yf,d.timeY);
    zrecon = interp1(t,zf,d.timeZ);
    
    figure()
    [hx, hy, hz] = plotConcentrations(d.timeX, d.xml, ...
        d.timeY, d.yml, d.timeZ, d.zml,'none');
    
    plotConcentrationLines(tr,xrepr,yrepr,tk,zk);
    
    dataRepresentation(i) = struct('data',d,...
        'xRepr',xrepr,'yRepr',yrepr,'zRepr',interp1(tk,zk,tr),...
        'alphaRepr',interp1(t,alpha{i},tr).','reprTime',tr.',...
        'xRecon',xrecon,'yRecon',yrecon,'zRecon',zrecon,...
        'alphaRecon',interp1(t,alpha{i},d.timeZ),'reconTime',d.timeZ,...
        'xExp',d.xml,'yExp',d.yml,'zExp',d.zml,...
        'alphaExp',alphaExp{i},'expTime',d.timeZ);
    
    save(savefilename);
    toc(tstart);
end
save(savefilename);
toc(tstart);
end

function plotConcentrationLines(tr,xr,yr,tk,zk)
subplot(1,3,1);
plotConcentrationLine(tr,xr);
xlabel('time [min]');
ylabel('concentration [mol/l]');
title('Alcohol-rich phase (X)');
subplot(1,3,2);
plotConcentrationLine(tr,yr);
xlabel('time [min]');
ylabel('concentration [mol/l]');
title('Oil-rich phase (Y)');
subplot(1,3,3);
plotConcentrationLine(tk,zk);
xlabel('time [min]');
ylabel('concentration [mol/l]');
title('Total mixture (Z)');
end

function plotConcentrationLine(t,f)
h = line(t,f);
set(h(1),'Color','b');
set(h(2),'Color','g');
set(h(3),'Color','r');
set(h(4),'Color','c');
set(h(5),'Color','m');
set(h(6),'Color','k');

set(h(1),'Marker','.');
set(h(2),'Marker','.');
set(h(3),'Marker','.');
set(h(4),'Marker','.');
set(h(5),'Marker','.');
set(h(6),'Marker','.');
end

function [xml2, yml2, mf] = SplitZToXY(xml,xTime, yml,yTime, zf, zT, alpha,time)
Vol = 1;
xVol = Vol*(1-alpha);
yVol = Vol*alpha;
xv = interp1(time,xVol,xTime);
yv = interp1(time,yVol,yTime);

for i =1:6
    zAtXTime(:,i) = interp1(zT,zf(:,i),xTime);
end

zm = zf.*repmat((xVol+yVol).',1,6);
zz = zAtXTime.*repmat((xv+yv),1,6);
ym = yml.*repmat(yv,1,6);
xm = xml.*repmat(xv,1,6);

dim = 4;
lBounds = -1e2*ones(dim,1);
uBounds = 1e2*ones(dim,1);
opts.MaxFunEvals = 1e6;
opts.Dim = dim;
p=2;
for (i = 1:6)
    [pars{i}, ~, ~, ~] = DeRandInfty(@GetSplittingError,...
        [], lBounds, uBounds, opts, xm(:,i), ym(:,i), zz(:,i), xTime,p);
    
    mf(:,i) = GetFunctionValue(pars{i},time,'gompertz');
    % mf must be between 0 and 1, if the regression method itself failed to
    % ensure it we need to force it (otherwise we'll have negative
    % concentrations)
    mf(:,i) = max(mf(:,i),0);
    mf(:,i) = min(mf(:,i),1);
    
    ym1(:,i) = mf(:,i).*zm(:,i);
    xm1(:,i) = (1-mf(:,i)).*zm(:,i);
end


xml2 = xm1./repmat(xVol.',1,6);
yml2 = ym1./repmat(yVol.',1,6);
% scaling = zf./(xf+yf);
% xf = xf.*scaling;
% yf = yf.*scaling;

if (norm(zf-repmat(alpha.',1,6).*yml2-(1-repmat(alpha.',1,6)).*xml2)>1e-6)
    warning('flash equation not met')
end


% figure()
% plot(yTime,y./zz,'o');
% hold on;
% legend('MeOH', 'GLY', 'TG', 'DG', 'MG', 'FAME')
% plot(zT,ym./(xm+ym),'.-');
% hold off;
end

function [err, x, y, mf] = GetSplittingError(pars,xm,ym,zm,time,p)
mf = GetFunctionValue(pars,time,'gompertz');
y = mf.*zm;
x = (1-mf).*zm;
err = norm(xm-x,p)+norm(ym-y,p);
end


function [xf, yf] = FitCurves(xml,xTime, yml,yTime, zml, zTime, xVol, yVol,time)
% MeOH, GLY, TG and FAME should be monotonous
for i = [1 2 3 6]
    [xf(:,i),yf(:,i)]=FitPair(xml(:,i),xTime,yml(:,i),yTime,zml(:,i),zTime, xVol, yVol, time,'gompertz');
end
% DG and MG need not to be monotonous
for i = [4 5]
    [xf(:,i),yf(:,i)]=FitPair(xml(:,i),xTime,yml(:,i),yTime,zml(:,i),zTime, xVol, yVol,time, 'bell');
end
end

function [alpha alphaExp] = ApproximateVolumes(xVol, xTime,yVol,yTime,time)
alphaExp = yVol./(xVol+yVol);
alpha = FitAlpha(alphaExp,xTime,time);

    function [v] = FitAlpha(vexp,texp,t)
        dim = 4;
        lBounds = -1e2*ones(dim,1);
        uBounds = 1e2*ones(dim,1);
        %         if (strcmp(monotonicity,'decr'))
        %             uBounds(2) = 0;
        %         elseif (strcmp(monotonicity,'incr'))
        %             lBounds(2) = 0;
        %         end
        %warning('Poor curve approximation');
        opts.MaxFunEvals = 1e6;
        opts.Dim = dim;
        err = @(pars,vexp,t) norm(vexp-GompertzFunction(pars, t),2);
        [pars, ~, ~, ~] = DeRandInfty(err,...
            [], lBounds, uBounds, opts, vexp,texp);
        v = GompertzFunction(pars, t);
    end

figure()
plot(xTime,alphaExp,'bo',time,alpha,'b');
legend('\alpha experimental', '\alpha approximation','Location','SE');
xlabel('time [min]');
ylabel('volume [ml]');
drawnow;
% hl = line(data{i}.timeX,1e2*data{i}.xVolume./(data{i}.xVolume+data{i}.yVolume));
% set(hl,'Marker','o');

% axis([0 70 0 1e2]);
% xlabel('time [min]');
% ylabel('Amount of oil phase [Vol. %]');
end

function [xf,yf]=FitPair(xml, xTime,yml,yTime,zml,zTime,xVol, yVol,time, type)
dim = 8;
lBounds = -10*ones(dim,1);
uBounds = 10*ones(dim,1);
opts.MaxFunEvals = 1e6;
opts.Dim = dim;
[pars, fMin, flag, output] = DeRandInfty(@CalculateFitError,...
    [], lBounds, uBounds, opts, xml, xTime, yml, yTime, zml, zTime, xVol, yVol, time, type);
[xf, yf] = GetFunctionValues(pars, time, type);
end

function v = GompertzFunction(pars, t)
a = pars(1);
b = pars(2);
c = -abs(log(abs(pars(3))));
d = -abs(log(abs(log(abs(pars(4))))));
v = a + b * exp(c*exp(d*t));
end

function v = BellFunction(pars,t)
a = pars(1);
b = pars(2);
c = log(abs(pars(3)));
d = log(abs(pars(4)));
v = a + b * exp(-((t-c).^2)./d);
end

function [xf, yf] = GetFunctionValues(pars, t, type)
xPars = pars(1:4);
yPars = pars(5:8);
xf = GetFunctionValue(xPars, t, type);
yf = GetFunctionValue(yPars, t, type);
end

function v = GetFunctionValue(pars, t, type)
if (strcmp(type,'gompertz'))
    v = GompertzFunction(pars,t);
elseif (strcmp(type,'bell'))
    v = BellFunction(pars,t);
else
    warning('CreateDataRepresentation:GetFunctionValues',...
        'Unsuported function type');
end
end

function err = CalculateFitError(pars, xml, xTime, yml, yTime, zml, zTime, xVol, yVol, time, type)
xf = GetFunctionValue(pars, xTime, type);
yf = GetFunctionValue(pars, yTime, type);
xv = interp1(time,xVol,zTime);
yv = interp1(time,yVol,zTime);

p = 0.5;
zv = xv+yv;
zf = xf.*(xv./zv)+yf.*(yv./zv);
error = [
    norm(xml-xf,p) % fit error for oil-rich phase
    norm(yml-yf,p) % fit error for alcohol-rich phase
    norm(zml-zf,p) % mass balance between x, y and z
    ];
penaltyForNegativeConc = sum(abs(xf(xf<0)))+sum(abs(yf(yf<0)));
err = norm(error,p)+10*norm(penaltyForNegativeConc,p);
end

function [xr, yr, tr] = TakeNewData(xf, yf, tf, N)
dx = diff(xf);
dy = diff(yf);
td = tf(1:(length(tf)-1))+0.5*mean(diff(tf));

for i = 1:length(td)
    delta(i) = norm(dx(i,:),2) + norm(dy(i,:),2);
end

cs = cumsum(delta);
cs = cs-cs(1);
empiricalCDF = cs./cs(length(cs));


toSample = linspace(min(empiricalCDF),max(empiricalCDF),N);

tr = NaN(1,N);
j = 1;
for i=1:N
    qi = toSample(i);
    while(qi<empiricalCDF(j) || qi>empiricalCDF(j+1))
        j=j+1;
    end
    cdPrev = empiricalCDF(j);
    cdFoll = empiricalCDF(j+1);
    tr(i) = td(j)+(td(j+1)-td(j))*(qi-cdPrev)./(cdFoll-cdPrev);
end
for i = 1:6
    xr(:,i) = interp1(tf,xf(:,i),tr);
    yr(:,i) = interp1(tf,yf(:,i),tr);
end
end






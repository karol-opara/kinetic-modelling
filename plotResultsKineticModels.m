function plotResultsKineticModels()
%plotExampleModel();
%plotZoomedMolels();
%plotTooFastReaction();

plotBenchmarkExample('Noureddini');
plotBenchmarkExample('Jansri');


%plotMotivatoryExample();
end

function plotBenchmarkExample(id)
% id = 'Jansri';
rndErr = 0.05;
sysErr = [1 0 0 0 0 0];



zc=[];
N=5;
for ii=1:N
    data{ii} = CreateBenchmarkProblem20130723(rndErr, sysErr, id);
    zc=[zc; data{ii}.zml; NaN(1,6)];
end

[tk, zk] = SolveKineticModel(data{1}.z0ml,data{1}.k,'batch');

ii=1;
lineStyle='-';
ls2 = ':';

if strcmp(id, 'Jansri')
    hu = subplot(2,2,2);
else
    hu = subplot(2,2,1);
end

hold on
ht1=plot(tk,zk(:,1),['b' lineStyle],tk,zk(:,2),['g' lineStyle],tk,zk(:,3),...
    ['r' lineStyle],tk,zk(:,4),['c' lineStyle],tk,zk(:,5),['m' lineStyle],...
    tk,zk(:,6),['k' lineStyle]);
set(ht1,'LineWidth',2);
zTime=repmat([data{ii}.timeZ NaN],1,N);
plot(zTime, zc(:,1), ['bd' ls2],zTime, zc(:,2), ['gx' ls2],zTime, zc(:,3), ['ro' ls2],...
    zTime, zc(:,4), ['c+' ls2],zTime, zc(:,5), ['m*' ls2],zTime, zc(:,6), ['ks' ls2]);
ax1 = [-0.5 20.5 -0.2 6];
axis(ax1)
xlabel('time [min]');
ylabel('concentration [mol/l]')
box('on');
if strcmp(id, 'Jansri')
    title('c)');
else
    title('a)');
end

if strcmp(id, 'Jansri')
    hd = subplot(2,2,4);
else
    hd = subplot(2,2,3);
end
hold on
hm=plot(tk,zk(:,1),['b' lineStyle],tk,zk(:,2),['g' lineStyle],tk,zk(:,3),...
    ['r' lineStyle],tk,zk(:,4),['c' lineStyle],tk,zk(:,5),['m' lineStyle],...
    tk,zk(:,6),['k' lineStyle]);
set(hm,'LineWidth',2);
hp = plot(zTime, zc(:,1), ['bd' ls2],zTime, zc(:,2), ['gx' ls2],zTime, zc(:,3), ['ro' ls2],...
    zTime, zc(:,4), ['c+' ls2],zTime, zc(:,5), ['m*' ls2],zTime, zc(:,6), ['ks' ls2]);
if strcmp(id, 'Jansri')
    ax2 = [-0.1 8 -0.02 3];
else
    ax2 = ax1;
    ax2(4) = 2.5;
end

axis(ax2)
xlabel('time [min]');
ylabel('concentration [mol/l]')
box('on');
if strcmp(id, 'Jansri')
    title('d)');
else
    title('b)');
end
if strcmp(id, 'Jansri')
    hleg = legend(hp,'MeOH', 'GLY', 'TG', 'DG', 'MG', 'FAME');
    set(hleg,'Orientation','Horizontal');
    pos = get(hleg,'Position');
    pos(2) = 0.001;
    pos(1) = (1-pos(3))/2;
    set(hleg,'Position',pos);
end

[x1, y1] = ds2nfu(hu, [ax1(1) ax2(2)],[ax1(3) ax1(3)]);
[x2, y2] = ds2nfu(hd, [ax2(1) ax2(2)],[ax2(4) ax2(4)]);
for i=1:numel(x1)
    h = annotation('line',[x1(i) x2(i)],[y1(i) y2(i)]);
    set(h,'LineStyle','--');
end

set(gcf, 'Position', [994   146   900   648]);


    function data = CreateBenchmarkProblem20130723(rndErr, sysErr, id)
        % first
        if (strcmp(id,'Oh'))
            k = [0.0500
                0.1100
                0.2150
                1.2280
                0.2420
                0.0070];
        elseif (strcmp(id,'Jansri'))
            k = [2.6
                0.248
                1.186
                0.227
                2.303
                0.022];
        elseif (strcmp(id,'Noureddini'))
            k = [0.0500
                0.1100
                0.2150
                1.2280
                0.2420
                0.0070];
        elseif (strcmp(id,'Klofutar'))
            k = [0.0443
                0.2334
                0.0645
                0.0699
                0.2681
                0.0047];
        else
            error('runBenchmarkingExperiment:WrongID','Wrong ID of the reaction rates');
        end
        
        conc = [6 0 1 0 0 0];
        
        % warning('Total mole per litre setting changed')
        % totalMolePerLitre = 10; % in simulations and initial submission
        totalMolePerLitre = 0.789359594 + 0.057795443 + 0.003075868 + 0.003075868 + 5.002450688; % in revision
        
        z0ml = totalMolePerLitre * conc/sum(conc);
        if (strcmp(id,'Oh'))
            timeZ = [2     4     6     8    10    12    20    30    40    50    60];
        elseif(strcmp(id,'Jansri'))
            timeZ = [0.5, 1, 3, 5, 7, 9, 12, 15, 18, 20 ];
        elseif (strcmp(id,'Noureddini'))
            timeZ = [1 2 3 4 5 6 8 10 15 20 30 50 60 90]; % Nouredini and Zhou (1997)
        end
        [t, z] = SolveKineticModel(z0ml,k);
        zKinetic = interp1(t,z,timeZ);
        
        rnd = rndErr* sqrt(pi/2) * randn(length(timeZ),length(z0ml));
        zml = zKinetic+zKinetic.*rnd+repmat(sysErr,length(timeZ),1);
        
        data = struct('k',k,'z0ml',z0ml,'timeZ',timeZ,'zml',zml,'zKinetic',zKinetic);
        
    end
end



function plotTooFastReaction()
ax3 = getAxis(0.1);
ax3(1:2) = [-0.2 5];

subplot(1,2,1)
hl1 = plotComparison('2','PoorDataNoRegularization');
axis(ax3);
title('a)')

subplot(1,2,2)
hl2 = plotComparison('2log','PoorDataNoRegularization');
axis(ax3);
title('b)');
end

function plotZoomedMolels()
figure()
ax1 = getAxis(10);
ax2 = getAxis(2.5);
ax3 = getAxis(0.15);
plotModelComparison('GoodDataNoRegularization', ax1, ax2, ax3);
figure()
plotModelComparison('GoodDataRegularization', ax1, ax2, ax3);

ax1 = getAxis(15);
ax2 = getAxis(0.8);
ax3 = getAxis(0.1);
figure()
plotModelComparison('PoorDataNoRegularization', ax1, ax2, ax3);
figure()
plotModelComparison('PoorDataRegularization', ax1, ax2, ax3);
end

function ax = getAxis(yMax)
ax = [-1 65 -yMax/30 yMax];
end

function plotExampleModel()
load('Results/save_LowConcComparison_2013-07-05_121424_reg01')
m = model{3,7};
ht = subplot(2,1,1);
plotKineticModelFit(m);
ax1 = [-0.7 65 -0.1 15];
axis(ax1);
title('a)')
legend('off')
hb = subplot(2,1,2);
[h,t,z,minus,hLegend] = plotKineticModelFit(m);
ax2 = [-0.7 65 -0.01 1.5];
axis(ax2);
title('b)')
set(gcf,'Position',[  9    49   467   919]);
set(hLegend, 'Orientation', 'Horizontal');
pos = get(hLegend,'Position');
pos(2) = 0.01;
pos(1) = (1-pos(3))/2;
set(hLegend,'Position',pos);

[x1, y1] = ds2nfu(ht, [ax1(1) ax1(2)],[ax1(3) ax2(3)]);
[x2, y2] = ds2nfu(hb, [ax2(1) ax2(2)],[ax2(4) ax2(4)]);

for i=1:numel(x1)
    h = annotation('line',[x1(i) x2(i)],[y1(i) y2(i)]);
    set(h,'LineStyle','--');
end

set(gcf, 'Position', [994   146   509   748])
end

function plotModelComparison(name,ax1,ax2,ax3,titleName)
set(gcf,'Position',[ 9    49   824   918]);
ht1 = subplot(3,2,1);
hl = plotComparison('2',name);
set(hl,'Visible','off');
axis(ax1);
title('a)')

hm1 = subplot(3,2,3);
hl = plotComparison('2',name);
set(hl,'Visible','off');
axis(ax2);
title('c)')
ch = get(gcf,'Children');
%set(ch(length(ch)-1),'Orientation','horizontal');

hb1 = subplot(3,2,5);
hl = plotComparison('2',name);
set(hl,'Visible','off');
axis(ax3);
title('e)')

ht2=subplot(3,2,2);
hl = plotComparison('2log',name);
set(hl,'Visible','off');
axis(ax1);
title('b)')

hm2=subplot(3,2,4);
hl = plotComparison('2log',name);
set(hl,'Visible','off');
axis(ax2);
title('d)')

hb2=subplot(3,2,6);
hl = plotComparison('2log',name);
set(hl,'Orientation','horizontal');
pos = get(hl,'Position');
pos(2) = 0.02;
pos(1) = (1-pos(3))/2;
set(hl,'Position',pos);
axis(ax3);
title('f)')

[xt1, yt1] = ds2nfu(ht2, [ax1(1) ax1(2)],[ax1(3) ax1(3)]);
[xt2, yt2] = ds2nfu(hm2, [ax2(1) ax2(2)],[ax2(4) ax2(4)]);
[xa1, ya1] = ds2nfu(hm2, [ax2(1) ax2(2)],[ax2(3) ax2(3)]);
[xa2, ya2] = ds2nfu(hb2, [ax3(1) ax3(2)],[ax3(4) ax3(4)]);
x1 = [xt1, NaN, xa1];
x2 = [xt2, NaN, xa2];
y1 = [yt1, NaN, ya1];
y2 = [yt2, NaN, ya2];
[xt1, yt1] = ds2nfu(ht1, [ax1(1) ax1(2)],[ax1(3) ax1(3)]);
[xt2, yt2] = ds2nfu(hm1, [ax2(1) ax2(2)],[ax2(4) ax2(4)]);
[xa1, ya1] = ds2nfu(hm1, [ax2(1) ax2(2)],[ax2(3) ax2(3)]);
[xa2, ya2] = ds2nfu(hb1, [ax3(1) ax3(2)],[ax3(4) ax3(4)]);
x1 = [x1, NaN, xt1, NaN, xa1];
x2 = [x2, NaN, xt2, NaN, xa2];
y1 = [y1, NaN, yt1, NaN, ya1];
y2 = [y2, NaN, yt2, NaN, ya2];
for i=1:numel(x1)
    if (~any(isnan([x1(i) x2(i) y1(i) y2(i)])))
        h = annotation('line',[x1(i) x2(i)],[y1(i) y2(i)]);
        set(h,'LineStyle','--');
    end
end
end

function hl =plotComparison(which,name)
if (strcmp(name,'GoodDataNoRegularization'))
    load('Results/save_LowConcComparison_2013-07-05_100552_reg0');
    m2 = model{1,2};
    m2log = model{3,2};
elseif (strcmp(name,'GoodDataRegularization'))
    load('Results/save_LowConcComparison_2013-07-05_183649_reg01');
    m2 = model{1,2};
    m2log = model{3,2};
elseif (strcmp(name,'PoorDataNoRegularization'))
    load('Results/save_LowConcComparison_2013-07-05_153144_reg0');
    m2 = model{1,4};
    m2log = model{4,4};
elseif (strcmp(name,'PoorDataRegularization'))
    load('Results/save_LowConcComparison_2013-07-05_183649_reg01');
    m2 = model{1,4};
    m2log = model{4,4};
end
if(strcmp(which,'2'))
    [~, ~,~, ~,hl]=plotKineticModelFit(m2);
elseif (strcmp(which,'2log'))
    [~, ~,~, ~,hl]=plotKineticModelFit(m2log);
end
end

function plotMotivatoryExample()
ax2 = getAxis(0.8);
subplot(1,2,1)
hl = plotComparison('2', 'PoorDataNoRegularization');
set(hl,'Visible','off');
axis(ax2);
title('a)')
gcapos = get(gca, 'Position');
shift = 0.1;
gcapos(2) = gcapos(2) + shift;
gcapos(4) = gcapos(4) - shift;
set(gca, 'Position', gcapos);
subplot(1,2,2)
hl = plotComparison('2log', 'PoorDataRegularization');
axis(ax2);
title('b)')
gcapos = get(gca, 'Position');
shift = 0.1;
gcapos(2) = gcapos(2) + shift;
gcapos(4) = gcapos(4) - shift;
set(gca, 'Position', gcapos);
set(hl,'Orientation','horizontal');
pos = get(hl,'Position');
pos(2) = 0.02;
pos(1) = (1-pos(3))/2;
set(hl,'Position',pos);
end
function plotResultsKineticModels()
%plotExampleModel()
%plotZoomedMolels();
plotTooFastReaction();
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
subplot(2,1,1)
plotKineticModelFit(m);
axis([0 65 0 15]);
title('a)')
subplot(2,1,2)
plotKineticModelFit(m);
axis([0 65 0 1.5]);
title('b)')
set(gcf,'Position',[  9    49   467   919]);
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
set(ch(length(ch)-1),'Orientation','horizontal');

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
    h = annotation('line',[x1(i) x2(i)],[y1(i) y2(i)]);
    set(h,'LineStyle','--');
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
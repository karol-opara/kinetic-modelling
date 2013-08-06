function [h,t,z,minus,hLegend]=plotKineticModelFit(zTime,zml,k,z0,type,breakAxis,lineStyle)
% PLOTKINETICMODELFIT plots the results of fitting kinetics model
% usage:
% model = EstimateKineticModel();
% plotKineticModelFit(model);
if (nargin == 1)
    model = zTime;
    zTime = model.data.timeZ;
    zml = model.data.zml;
    k = model.k;
    z0 = model.z0opt;
    type = model.type;
    breakAxis = false;
end
if (nargin < 6)
    breakAxis = false;
    lineStyle = '-';
end

[t, z] = SolveKineticModel(z0,k,type);

h=plot(zTime, zml(:,1), 'bd',zTime, zml(:,2), 'gx',zTime, zml(:,3), 'ro',...
    zTime, zml(:,4), 'c+',zTime, zml(:,5), 'm*',zTime, zml(:,6), 'ks',...
    t,z(:,1),['b' lineStyle],t,z(:,2),['g' lineStyle],t,z(:,3),...
    ['r' lineStyle],t,z(:,4),['c' lineStyle],t,z(:,5),['m' lineStyle],...
    t,z(:,6),['k' lineStyle]);
for i=1:length(h)
    set(h,'LineWidth',1);
end
xlabel('time [min]');
%legend('MeOH', 'GLY', 'TG', 'DG', 'MG','FAME','Orientation','Horizontal')
hLegend = legend('MeOH', 'GLY', 'TG', 'DG', 'MG','FAME');
xlabel('time [min]');
ylabel('concentration [mol/l]')
axis([-1 70.2 -0.2 Inf])
set(gca,'TickDir','out');
drawnow;

ystart = ceil(max(max([zml(:,2:6); z(:,2:6)])));
ystop = floor(min(min([zml(:,1); z(:,1)])));
if (breakAxis)
    minus=MyBreakAxis(gca,h,ystart,ystop);
else
    minus = 0;
end
end

function minus = MyBreakAxis(gca, hl, ystart, ystop)
width = 1;
minus = (ystop - ystart - width);
for i =1:length(hl)
    y = get(hl(i),'YData');
    y(y>ystart) = y(y>ystart) - minus;
    set(hl(i),'YData',y);
end
xtick=get(gca,'XTick');
t1=text(xtick(1), ystart+width/2,'//','fontsize',15);
t2=text(xtick(max(length(xtick))),ystart+width/2,'//','fontsize',15);
set(t1,'rotation',270);
set(t2,'rotation',270);

% remap tick marks, and 'erase' them in the gap
ytick=get(gca,'yTick');
dtick=ytick(2)-ytick(1);
gap=floor(width/dtick);
last=max(ytick(ytick<=ystart));          % last tick mark in LH dataset
next=min(ytick(ytick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
offset = width;
%offset=size(y2(y2>last&y2<next),2)*(y(2)-y(1));

for i=1:sum(ytick>(last+gap))
    ytick(find(ytick==last)+i+gap)=ystop+offset+dtick*(i-1);
end

lastI = 1;
for i=1:length(ytick)
    if ytick(i)>last&ytick(i)<next
        yticklabel{i}=sprintf('%d',[]);
        lastI = i;
    else
        yticklabel{i}=num2str(ytick(i));
    end
end;
yticklabel{lastI} = num2str(ytick(lastI+1)-(ytick(lastI+2)-(ytick(lastI+1))));
set(gca,'yticklabel',yticklabel);

end
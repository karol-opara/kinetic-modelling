function [hx, hy, hz] = plotConcentrations(timeX, x, timeY, y, timeZ, z,lineStyle)
if (nargin < 7)
    lineStyle = '-';
end
if (strcmp(lineStyle,'none'))
    lineStyle = '';
end
hx = subplot(1,3,1);
plotConcentration(timeX, x, lineStyle);
hy = subplot(1,3,2);
plotConcentration(timeY, y, lineStyle);
hz = subplot(1,3,3);
plotConcentration(timeZ, z, lineStyle);
legend('MeOH', 'GLY', 'TG', 'DG', 'MG', 'FAME');
end

function [h] = plotConcentration(time, concentration, lineStyle)
h = plot(time, concentration(:,1), ['bo' lineStyle],...
    time, concentration(:,2), ['go' lineStyle],...
    time, concentration(:,3), ['ro' lineStyle],...
    time, concentration(:,4), ['co' lineStyle],...
    time, concentration(:,5), ['mo' lineStyle],...
    time, concentration(:,6), ['ko' lineStyle]);
marginY = 0;
marginX = 0;
yl = ylim;
if (yl(2)>30)
    ylim2 = max(yl(2),100);
else
    ylim2 = 30;
end
set(gca,'YGrid','on','YMinorGrid','on')
axis([-marginX, 70+marginX, -marginY, ylim2+marginY]);
end
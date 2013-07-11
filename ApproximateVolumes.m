function [xVol, yVol] = ApproximateVolumes(xVol,yVol,time)
load Results\saveExperimentalData20120807

for i=1:length(data)%[3 5 6 7]
    hl = line(data{i}.timeX,1e2*data{i}.xVolume./(data{i}.xVolume+data{i}.yVolume));
    set(hl,'Marker','o');
end
axis([0 70 0 1e2]);
xlabel('time [min]');
ylabel('Amount of oil phase [Vol. %]');
end
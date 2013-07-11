function processKineticsLpComparison()
load Results/saveKineticsModels_DifferentLpNormsWithLongOptimizationTime.mat
close all
[lenData, lenLp] = size(models);
lineStyles = {'', ':', '--', '-'};
for i = 1:lenData
    no = figure();
    set(no,'Position',[100 100 550 850]);
    m= models{i,1};
    [h,t,z,minus] = plotKineticModelFit(models{i,2}.data.timeZ, ...
        models{i,2}.data.zml, models{i,2}.k, models{i,2}.z0opt,...
        true,':');
    
    [t1,z1] = SolveKineticModel(models{i,3}.z0opt, models{i,3}.k);
    h1 = plotKinetics(t1,z1,minus,'--');
        
    [t2,z2] = SolveKineticModel(models{i,4}.z0opt, models{i,4}.k);
     h2 = plotKinetics(t2,z2,minus,'-');
     
     legend([h(1:6).' h(12) h1(6) h2(6)],...
         'MeOH', 'GLY', 'TG', 'DG', 'MG','FAME',...
         'L0.5 fit', 'L1 fit', 'L2 fit');
    legend('boxoff');
    
    for j = 1:lenLp
        m= models{i,j};

        K{i,j} = struct('k',m.k,'kf',sum(m.k([1 3 5])), 'kb',sum(m.k([2 4 6])),...
            'id',['ex ' num2str(i) ' L' num2str(p(j))]);
    end
end
%save saveExperimentalK_DifferentLpNormsWithLongOptimizationTime
end

function h = plotKinetics(t,z, minus, lineStyle)
hold on
h = plot(t,z(:,1)-minus,['b' lineStyle],t,z(:,2),['g' lineStyle],t,z(:,3),...
    ['r' lineStyle],t,z(:,4),['c' lineStyle],t,z(:,5),['m' lineStyle],...
    t,z(:,6),['k' lineStyle]);
hold off
end
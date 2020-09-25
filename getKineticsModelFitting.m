function getKineticsModelFitting
load Results\RegNonReg\saveKineticsModels_p2_fitReg

plotKineticModelCurves(data, models);
%compareFameConcentration(data, models);

end

function plotKineticModelCurves(data, models)
for i = 1:length(models)
    figure()
    plotKineticModelFit(data{i}.timeZ,data{i}.zml, models{i}.k, data{i}.z0ml);
    title(getExperimentDescription(i,data));
    axis([0 70 0 15]);
    set(gca,'YMinorGrid','on', 'YGrid','on')
end
end

function compareFameConcentration(data, models)
axis();
experimentsToCompare = [1 2 3 4 5 6 7];
colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
legendEntries = {};
for i = experimentsToCompare
    [t, z] = SolveKineticModel(data{i}.z0ml, models{i}.k);
    hl = line(t,z(:,6));
    set(hl,'Color',colors{i});
    hd = line(data{i}.timeZ,data{i}.zml(:,6));
    set(hd,'Color',colors{i},'Marker','o','MarkerFaceColor',colors{i},'LineStyle','none');
    legendHandles(i) = hl;
    legendEntries{i}= getExperimentDescription(i,data);
end
legend(legendHandles,legendEntries);
xlabel('time [min]');
ylabel('biodiesel concentration [mol/l]');
set(gca,'YMinorGrid','on', 'YGrid','on')
end

function description = getExperimentDescription(i,data)
description = ['exp. ' num2str(i)...
    ', NaOH = ' num2str(data{i}.NaOH)...
    ', TG:MeOH = ' num2str(data{i}.TGMeOHRatio)];
end
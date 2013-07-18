function processLandscapes()
processLandscapes2013();
%processLandscapes2012(); % landscapes from YSSP report
end

function processLandscapes2013()
poorDataFit={'Results/save_2013-07-11_163658_KineticsModelLandscape_condition_5_SRI_nonregularized_L2_QNaN',...
    'Results/save_2013-07-12_032246_KineticsModelLandscape_condition_5_SRI_regularized01_L2_QNaN',...
    'Results/save_2013-07-13_124040_KineticsModelLandscape_condition_5_SRI_nonregularized_L2_Qlog',...
    'Results/save_2013-07-14_000853_KineticsModelLandscape_condition_5_SRI_regularized01_L2_Qlog'};
goodDataFit={'Results/save_2013-07-12_133850_KineticsModelLandscape_condition_7_SRI_nonregularized_L2_QNaN',...
    'Results/save_2013-07-13_005848_KineticsModelLandscape_condition_7_SRI_regularized01_L2_QNaN',...
    'Results/save_2013-07-14_115358_KineticsModelLandscape_condition_7_SRI_nonregularized_L2_Qlog',...
    'Results/save_2013-07-14_234202_KineticsModelLandscape_condition_7_SRI_regularized01_L2_Qlog'};
l1data = {'Results/save_2013-07-15_113034_KineticsModelLandscape_condition_5_SRI_nonregularized_L1_Qlog',...
    'Results/save_2013-07-15_224208_KineticsModelLandscape_condition_5_SRI_regularized01_L1_Qlog',...
    'Results/save_2013-07-16_101005_KineticsModelLandscape_condition_7_SRI_nonregularized_L1_Qlog',...
    'Results/save_2013-07-16_220230_KineticsModelLandscape_condition_7_SRI_regularized01_L1_Qlog'};
    

close all
titles = {' L2 NaN',' regularized L2 NaN',' L2 log', ' regularized L2 log'};
plotFourLandscapes(poorDataFit,titles);
plotFourLandscapes(goodDataFit,titles);
titles = {' poor data L1 log', ' poor data regularizes L1 log', ' good data L1 log', ' good data regularized L1 log'};
plotFourLandscapes(l1data,titles);
end

function plotFourLandscapes(data,titles)
nr = figure();
subplot(2,2,1);
processLandscape(data{1},'contourf');
title(['a)' titles{1}]);
subplot(2,2,2);
processLandscape(data{2},'contourf');
title(['b)' titles{2}]);
subplot(2,2,3);
processLandscape(data{3},'contourf');
title(['c)' titles{3}]);
subplot(2,2,4);
processLandscape(data{4},'contourf');
title(['d)' titles{4}]);

set(nr,'Position',[100,100, 1000, 700']);
end

function processLandscapes2012()
pathsOrigWithoutGLY_Lp_Lp={'Results/save_2012-08-08_153735_KineticsModelLandscape_condition_5_withGLY_Lp025'...
    'Results/save_2012-08-08_153749_KineticsModelLandscape_condition_5_withGLY_Lp05'...
    'Results/save_2012-08-08_153759_KineticsModelLandscape_condition_5_withGLY_Lp1'...
    'Results/save_2012-08-08_153814_KineticsModelLandscape_condition_5_withGLY_Lp2'...
    'Results/save_2012-08-08_153839_KineticsModelLandscape_condition_7_withGLY_Lp025'...
    'Results/save_2012-08-08_153854_KineticsModelLandscape_condition_7_withGLY_Lp05'...
    'Results/save_2012-08-08_153942_KineticsModelLandscape_condition_7_withGLY_Lp1'...
    'Results/save_2012-08-08_153952_KineticsModelLandscape_condition_7_withGLY_Lp2'...
    };

pathWithGLY_Lp_L2 ={'Results\save_2012-08-14_150852_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L025',...
    'Results\save_2012-08-14_082006_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L05',...
    'Results\save_2012-08-14_082004_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L1',...
    'Results\save_2012-08-14_082000_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L2',...
    'Results\save_2012-08-14_150900_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L025',...
    'Results\save_2012-08-14_081954_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L05',...
    'Results\save_2012-08-14_081118_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L1',...
    'Results\save_2012-08-14_081115_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L2'};

close all

figure(1)
for i = 1:4
    subplot(2,2,i);
    processLandscape(pathWithGLY_Lp_L2{i},'contourf');
end
figure(2)
for i=1:4
    subplot(2,2,i);
    processLandscape(pathWithGLY_Lp_L2{i},'plot');
end

figure(3)
for i = 5:8
    subplot(2,2,i-4);
    processLandscape(pathWithGLY_Lp_L2{i},'contourf');
end
figure(4)
for i = 5:8
    subplot(2,2,i-4);
    processLandscape(pathWithGLY_Lp_L2{i},'plot');
end

for i = 1:4
    set(i,'Position',[100,100, 1000, 700']);
end

end
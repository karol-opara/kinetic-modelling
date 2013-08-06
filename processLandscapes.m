function processLandscapes()
processLandscapes2013();
%processLandscapes2012(); % landscapes from YSSP report
end

function processLandscapes2013()
% poorDataFit={'Results/save_2013-07-11_163658_KineticsModelLandscape_condition_5_SRI_nonregularized_L2_QNaN',...
%     'Results/save_2013-07-12_032246_KineticsModelLandscape_condition_5_SRI_regularized01_L2_QNaN',...
%     'Results/save_2013-07-13_124040_KineticsModelLandscape_condition_5_SRI_nonregularized_L2_Qlog',...
%     'Results/save_2013-07-14_000853_KineticsModelLandscape_condition_5_SRI_regularized01_L2_Qlog'};
% goodDataFit={'Results/save_2013-07-12_133850_KineticsModelLandscape_condition_7_SRI_nonregularized_L2_QNaN',...
%     'Results/save_2013-07-13_005848_KineticsModelLandscape_condition_7_SRI_regularized01_L2_QNaN',...
%     'Results/save_2013-07-14_115358_KineticsModelLandscape_condition_7_SRI_nonregularized_L2_Qlog',...
%     'Results/save_2013-07-14_234202_KineticsModelLandscape_condition_7_SRI_regularized01_L2_Qlog'};
% l1data = {'Results/save_2013-07-15_113034_KineticsModelLandscape_condition_5_SRI_nonregularized_L1_Qlog',...
%     'Results/save_2013-07-15_224208_KineticsModelLandscape_condition_5_SRI_regularized01_L1_Qlog',...
%     'Results/save_2013-07-16_101005_KineticsModelLandscape_condition_7_SRI_nonregularized_L1_Qlog',...
%     'Results/save_2013-07-16_220230_KineticsModelLandscape_condition_7_SRI_regularized01_L1_Qlog'};
%
% l2nonregData = {'Results/save_2013-07-23_095951_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_QNaN',...
%     'Results/save_2013-07-23_232017_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_QNaN',...
%     'Results/save_2013-07-24_132350_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_Qlog',...
%     'Results/save_2013-07-25_033347_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_Qlog'};
%
% l1nonregData = {'Results/save_2013-07-25_174438_KineticsModelLandscape_condition_5_OPP_nonregularized_L1_Qlog',...
%     'Results/save_2013-07-26_080519_KineticsModelLandscape_condition_7_OPP_nonregularized_L1_Qlog',...
%     'Results/save_2013-07-26_223843_KineticsModelLandscape_condition_5_OPP_nonregularized_Llog_QNaN',...
%     'Results/save_2013-07-27_115222_KineticsModelLandscape_condition_7_OPP_nonregularized_Llog_QNaN'};

L2QNaNData = {'Results/save_2013-07-23_095951_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_QNaN',... % 5 non
    'Results/save_2013-07-29_100926_KineticsModelLandscape_condition_5_OPP_regularized_L2_QNaN',... % 5 reg
    'Results/save_2013-07-23_232017_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_QNaN',... % 7 non
    'Results/save_2013-07-29_222219_KineticsModelLandscape_condition_7_OPP_regularized_L2_QNaN'}; % 7 reg

L2QLogData = {'Results/save_2013-07-24_132350_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_Qlog',... % 5 non
    'Results/save_2013-07-30_110913_KineticsModelLandscape_condition_5_OPP_regularized_L2_Qlog',... % 5 reg
    'Results/save_2013-07-25_033347_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_Qlog',... % 7 non
    'Results/save_2013-07-31_001330_KineticsModelLandscape_condition_7_OPP_regularized_L2_Qlog'}; % 7 reg

L1QLogData = {'Results/save_2013-07-25_174438_KineticsModelLandscape_condition_5_OPP_nonregularized_L1_Qlog',... % 5 non
    'Results/save_2013-07-31_130852_KineticsModelLandscape_condition_5_OPP_regularized_L1_Qlog',... % 5 reg
    'Results/save_2013-07-26_080519_KineticsModelLandscape_condition_7_OPP_nonregularized_L1_Qlog',... % 7 non
    'Results/save_2013-08-01_021127_KineticsModelLandscape_condition_7_OPP_regularized_L1_Qlog'}; % 7 reg

LLogQNaNData = {'Results/save_2013-07-26_223843_KineticsModelLandscape_condition_5_OPP_nonregularized_Llog_QNaN',... % 5 non
    'Results/save_2013-08-01_151919_KineticsModelLandscape_condition_5_OPP_regularized_Llog_QNaN',... % 5 reg
    'Results/save_2013-07-27_115222_KineticsModelLandscape_condition_7_OPP_nonregularized_Llog_QNaN',... % 7 non
    'Results/save_2013-08-02_041336_KineticsModelLandscape_condition_7_OPP_regularized_Llog_QNaN'}; % 7 reg


close all
titles = {'','','',''};%{' poor', ' poor regularized', ' good', ' good regularized'};
plotFourLandscapes(L2QNaNData,titles,'contourf');
plotFourLandscapes(L2QLogData,titles,'contourf');
%plotFourLandscapes(L1QLogData,titles,'contourf');
%plotFourLandscapes(LLogQNaNData,titles,'contourf');

plotFourLandscapes(L2QNaNData,titles,'plot');
plotFourLandscapes(L2QLogData,titles,'plot');
%plotFourLandscapes(L1QLogData,titles,'plot');
%plotFourLandscapes(LLogQNaNData,titles,'plot');

% titles = {' L2 NaN',' regularized L2 NaN',' L2 log', ' regularized L2 log'};
% plotFourLandscapes(poorDataFit,titles);
% plotFourLandscapes(goodDataFit,titles);
% titles = {' poor data L1 log', ' poor data regularizes L1 log', ' good data L1 log', ' good data regularized L1 log'};
% plotFourLandscapes(l1data,titles);
%
% titles={' poor L2 QNaN', ' good L2 QNaN', ' poor L2 QLog', ' good L2 QLog'};
% plotFourLandscapes(l2nonregData,titles);
%
% titles={' poor L1 Qlog', ' good L1 Qlog', ' poor Llog QNaN', ' good Llog QNan'};
% plotFourLandscapes(l1nonregData,titles);
end

function plotFourLandscapes(data,titles, type)
nr = figure();
subplot(2,2,1);
hl = processLandscape(data{1},type);
if(isnan(hl)==false), set(hl,'Visible','off'); end
title(['a)' titles{1}]);
subplot(2,2,2);
hl = processLandscape(data{2},type);
if(isnan(hl)==false), set(hl,'Visible','off'); end
title(['b)' titles{2}]);
subplot(2,2,3);
hl = processLandscape(data{3},type);
if(isnan(hl)==false), set(hl,'Visible','off'); end
title(['c)' titles{3}]);
subplot(2,2,4);
hl = processLandscape(data{4},type);
if(isnan(hl)==false)
    set(hl,'Orientation','horizontal');
    pos = get(hl,'Position');
    pos(2) = 0.005;
    pos(1) = (1-pos(3))/2;
    set(hl,'Position',pos);
end
title(['d)' titles{4}]);

set(nr,'Position',[100,100, 700, 700]);
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
function processLandscapes()
processLandscapes2013();
end

function processLandscapes2013()

L2QNaNData = {'../data/save_2013-07-23_095951_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_QNaN',... % 5 non
    '../data/save_2013-07-29_100926_KineticsModelLandscape_condition_5_OPP_regularized_L2_QNaN',... % 5 reg
    '../data/save_2013-07-23_232017_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_QNaN',... % 7 non
    '../data/save_2013-07-29_222219_KineticsModelLandscape_condition_7_OPP_regularized_L2_QNaN'}; % 7 reg

L2QLogData = {'../data/save_2013-07-24_132350_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_Qlog',... % 5 non
    '../data/save_2013-07-30_110913_KineticsModelLandscape_condition_5_OPP_regularized_L2_Qlog',... % 5 reg
    '../data/save_2013-07-25_033347_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_Qlog',... % 7 non
    '../data/save_2013-07-31_001330_KineticsModelLandscape_condition_7_OPP_regularized_L2_Qlog'}; % 7 reg

L1QLogData = {'../data/save_2013-07-25_174438_KineticsModelLandscape_condition_5_OPP_nonregularized_L1_Qlog',... % 5 non
    '../data/save_2013-07-31_130852_KineticsModelLandscape_condition_5_OPP_regularized_L1_Qlog',... % 5 reg
    '../data/save_2013-07-26_080519_KineticsModelLandscape_condition_7_OPP_nonregularized_L1_Qlog',... % 7 non
    '../data/save_2013-08-01_021127_KineticsModelLandscape_condition_7_OPP_regularized_L1_Qlog'}; % 7 reg

LLogQNaNData = {'../data/save_2013-07-26_223843_KineticsModelLandscape_condition_5_OPP_nonregularized_Llog_QNaN',... % 5 non
    '../data/save_2013-08-01_151919_KineticsModelLandscape_condition_5_OPP_regularized_Llog_QNaN',... % 5 reg
    '../data/save_2013-07-27_115222_KineticsModelLandscape_condition_7_OPP_nonregularized_Llog_QNaN',... % 7 non
    '../data/save_2013-08-02_041336_KineticsModelLandscape_condition_7_OPP_regularized_Llog_QNaN'}; % 7 reg


close all
titles = {'','','',''};%{' poor', ' poor regularized', ' good', ' good regularized'};
p1 = plotFourLandscapes(L2QNaNData,titles,'contourf');
p2 = plotFourLandscapes(L2QLogData,titles,'contourf');

p3 = plotFourLandscapes(L2QNaNData,titles,'plot');
p4 = plotFourLandscapes(L2QLogData,titles,'plot');

saveas(p1, '../results/Fig_6_landscapes_square.png')
saveas(p2, '../results/Fig_8_landscapes_logsquare.png')

saveas(p3, '../results/Fig_7_best_trajectories_square.png')
saveas(p4, '../results/Fig_9_best_trajectories_logsquare.png')
end

function nr = plotFourLandscapes(data,titles, type)
nr = figure();
subplot(2,2,1);
hl = processLandscape(data{1},type);
if(isnan(double(hl))==false), set(hl,'Visible','off'); end
title(['a)' titles{1}]);
subplot(2,2,2);
hl = processLandscape(data{2},type);
if(isnan(double(hl))==false), set(hl,'Visible','off'); end
title(['b)' titles{2}]);
subplot(2,2,3);
hl = processLandscape(data{3},type);
if(isnan(double(hl))==false), set(hl,'Visible','off'); end
title(['c)' titles{3}]);
subplot(2,2,4);
hl = processLandscape(data{4},type, true);
%if(isnan(double(hl))==false)
    set(hl,'Orientation','horizontal');
    pos = get(hl,'Position');
    pos(2) = 0.001;
    pos(1) = (1-pos(3))/2;
    set(hl,'Position',pos);
%end
title(['d)' titles{4}]);

set(nr,'Position',[100,100, 700, 700]);
end

function processLandscapes2012()
pathsOrigWithoutGLY_Lp_Lp={'../data/save_2012-08-08_153735_KineticsModelLandscape_condition_5_withGLY_Lp025'...
    '../data/save_2012-08-08_153749_KineticsModelLandscape_condition_5_withGLY_Lp05'...
    '../data/save_2012-08-08_153759_KineticsModelLandscape_condition_5_withGLY_Lp1'...
    '../data/save_2012-08-08_153814_KineticsModelLandscape_condition_5_withGLY_Lp2'...
    '../data/save_2012-08-08_153839_KineticsModelLandscape_condition_7_withGLY_Lp025'...
    '../data/save_2012-08-08_153854_KineticsModelLandscape_condition_7_withGLY_Lp05'...
    '../data/save_2012-08-08_153942_KineticsModelLandscape_condition_7_withGLY_Lp1'...
    '../data/save_2012-08-08_153952_KineticsModelLandscape_condition_7_withGLY_Lp2'...
    };

pathWithGLY_Lp_L2 ={'../data/save_2012-08-14_150852_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L025',...
    '../data/save_2012-08-14_082006_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L05',...
    '../data/save_2012-08-14_082004_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L1',...
    '../data/save_2012-08-14_082000_KineticsModelLandscape_condition_5_AllCompsUnnormalizedSecondAggrL2_L2',...
    '../data/save_2012-08-14_150900_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L025',...
    '../data/save_2012-08-14_081954_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L05',...
    '../data/save_2012-08-14_081118_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L1',...
    '../data/save_2012-08-14_081115_KineticsModelLandscape_condition_7_AllCompsUnnormalizedSecondAggrL2_L2'};

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
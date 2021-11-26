function prvKineticLandscapesRegularized

l2nonregData = {'Results/save_2013-07-23_095951_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_QNaN',...
    'Results/save_2013-07-23_232017_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_QNaN',...
    'Results/save_2013-07-24_132350_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_Qlog',...
    'Results/save_2013-07-25_033347_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_Qlog'};

l1nonregData = {'Results/save_2013-07-25_174438_KineticsModelLandscape_condition_5_OPP_nonregularized_L1_Qlog',...
    'Results/save_2013-07-26_080519_KineticsModelLandscape_condition_7_OPP_nonregularized_L1_Qlog',...
    'Results/save_2013-07-26_223843_KineticsModelLandscape_condition_5_OPP_nonregularized_Llog_QNaN',...
    'Results/save_2013-07-27_115222_KineticsModelLandscape_condition_7_OPP_nonregularized_Llog_QNaN'};

relativeLambda = 0.1;
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-23_095951_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_QNaN');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-23_232017_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_QNaN');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-24_132350_KineticsModelLandscape_condition_5_OPP_nonregularized_L2_Qlog');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-25_033347_KineticsModelLandscape_condition_7_OPP_nonregularized_L2_Qlog');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-25_174438_KineticsModelLandscape_condition_5_OPP_nonregularized_L1_Qlog');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-26_080519_KineticsModelLandscape_condition_7_OPP_nonregularized_L1_Qlog');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-26_223843_KineticsModelLandscape_condition_5_OPP_nonregularized_Llog_QNaN');
runKineticRegularizedModel('OPP_regularized', relativeLambda, 'save_2013-07-27_115222_KineticsModelLandscape_condition_7_OPP_nonregularized_Llog_QNaN');
end

function runKineticRegularizedModel(name,relativeLambda, nonregularizedFilename)
[nonregularizedErrors, p, q, condition] = loadNonregularizedErrors(nonregularizedFilename);
runKinteticsModelLandscape(name, condition, p, q, relativeLambda, nonregularizedErrors)
end

function [err, p, q, condition] = loadNonregularizedErrors(nonregularizedFilename)
load(['Results/' nonregularizedFilename]);
[r, c] = size(Output);
err = NaN(r,c);
for i =1:r
    for j = 1:c
        err(i,j) = Output{i,j}.solutions.bestever.f;
    end
end
end
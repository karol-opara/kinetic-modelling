function processLowConcComparison
%processTwo();
processMany();
end

function processMany()
%load('Results/save_LowConcComparison_2013-07-05_100552_reg0'); % data{2} is as predicted
%load('Results/save_LowConcComparison_2013-07-05_111556_reg001');
%load('Results/save_LowConcComparison_2013-07-05_121424_reg01');

% all eight experiments
load('Results/save_LowConcComparison_2013-07-05_153144_reg0');
load('Results/save_LowConcComparison_2013-07-05_171134_reg001');
load('Results/save_LowConcComparison_2013-07-05_183649_reg01');
%load('Results/save_LowConcComparison_2013-07-05_193059_reg1');

for j = 4%1:8
    for i =[1 4]%1:4
        plotConcentrations3(model{i,j},['p = ' num2str(p{i}) ', q = ' num2str(q{i})]);
    end
end
end

function plotConcentrations3(model, name)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plotKineticModelFit(model);
title(name);
subplot(1,3,2)
plotKineticModelFit(model);
title(name);
axis([-1 70 -Inf 2.3]);
subplot(1,3,3)
plotKineticModelFit(model);
title(name);
axis([-1 70 -Inf 0.3]);
end

function processTwo()
load('Results/saveLowConcComparison') % with cmaes and 1e4 FEs
plotConcentrationsOKERR(m2_ok, m2_err,'m2');
plotConcentrationsOKERR(m1_ok, m1_err,'m1');
plotConcentrationsOKERR(m05_ok, m05_err,'m05');
plotConcentrationsOKERR(mlog_ok, mlog_err,'mlog');
plotConcentrationsOKERR(mlog_1_ok, mlog_1_err,'mlog 1');
plotConcentrationsOKERR(mlog_2_ok, mlog_2_err,'mlog 2');
end

function plotConcentrationsOKERR(m_ok, m_err, name)
figure()
subplot(2,3,1)
plotKineticModelFit(m_ok);
title([name ' ok']);
subplot(2,3,4)
plotKineticModelFit(m_err);
title([name ' err']);
subplot(2,3,2)
plotKineticModelFit(m_ok);
title([name ' ok']);
axis([-Inf 70 -Inf 2.3]);
subplot(2,3,5)
plotKineticModelFit(m_err);
title([name ' err']);
axis([-Inf 70 -Inf 1.7]);
subplot(2,3,3)
plotKineticModelFit(m_ok);
title([name ' ok']);
axis([-Inf 70 -Inf 0.2]);
subplot(2,3,6)
plotKineticModelFit(m_err);
title([name ' err']);
axis([-Inf 70 -Inf 0.2]);
end
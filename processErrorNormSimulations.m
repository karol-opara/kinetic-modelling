function good = processErrorNormSimulations

good = processRunUniqunessExperimentRepetetiveFits()
% processUniqunessExperimentImportanceSampling()

end

function good = processRunUniqunessExperimentRepetetiveFits()
% 10 reps
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-07-02_175030_IIASA_KO_random_05_sys_0_regularized_001_Nouredini_10rep'); % OK, 10% accuracy
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-07-03_063951_IIASA_KO_random_05_sys_0_regularized_001_Jansri_10rep'); % poor results but only outer log may help
% load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-07-03_102426_IIASA_PPO_random_05_sys_0_regularized_0001_Jansri_10rep'); % too poor regularization, no good results
% load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-07-03_211010_IIASA_KO_random_05_sys_0_regularized_01_Jansri_10rep');

% 25 reps
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-06-17_172155_IIASA_random_15_sys_0_regularized_001');
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-06-17_151119_SRI_random_15_sys_3_regularized_001')


% 100 reps
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-06-20_173831_SRI_random_15_sys_0_regularized_001')
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-06-24_230921_SRI_random_15_sys_3_regularized_001')
%load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-06-26_174839_IIASA_KO_random_05_sys_0_regularized_001')


% 100 reps with log regularization
load('Results/save_ErrorNormSimulations_RepetetiveFits_2013-07-09_172446_IIASA_PPO_random_05_sys_0_regularized_01');

lRelativeAccuracy = 0.5;
uRelativeAccuracy = 1.5;

if(false == exist('plen'))
   plen = len;
   qlen = len;
   pnorms = norms;
   qnorms = norms;
else
   rep = N;
end

kVal = cell(3,3);
for i = 1:plen
    for j = 1:qlen
        mtr = NaN(rep,6);
        for (r = 1:rep)
            mtr(r,:) = models{r,i,j}.k.';
        end
        kVal{i,j}= mtr;
    end
end
good=NaN;
for i = 1:plen
    for j = 1:qlen
        kf = sum(kVal{i,j}(:,[1 3 5]).').';
        kb = sum(kVal{i,j}(:,[2 4 6]).').';
        subplot(plen,qlen,qlen*(i-1)+j)
        h = plot(kf,kb,'go',sum(data.k([1,3,5])),sum(data.k([2,4,6])),'*');
        set(gca,'xScale','log','yScale','log');
        %axis([1e-3 100 1e-3 200]);
        axis([1e-2 1e3 1e-2 1e3]);
        title(['p = ' num2str(pnorms{i}) ', q = ' num2str(qnorms{j})]);
        xlabel('k_f = k_1 + k_3 +k_5')
        ylabel('k_b = k_2 + k_4 +k_6')
    
        [r, ~] = size(kVal{i,j});
        kc = kVal{i,j}(all(kVal{i,j}>repmat(data.k.'*lRelativeAccuracy,r,1),2) & ...
            all(kVal{i,j}<repmat(data.k.'*uRelativeAccuracy,r,1),2),:);
        hl = line(sum(kc(:,[1 3 5]).').',sum(kc(:,[2 4 6]).').');
        set(hl,'color','r','marker','o','linestyle','none');
        good(i,j) = round(100*numel(kc)/numel(kVal{i,j}));
        nvar(i,j) = norm(kVal{i,j});
        title(['p = ' num2str(pnorms{i}) ', q = ' num2str(qnorms{j}) ', good = ' ...
            num2str(good(i,j)) '%']);
    end
end
good
%round(nvar)
%figure()
%bar3(nvar)
%bar3(good)

end

function processUniqunessExperimentImportanceSampling()
names={'OnlyRandom', 'RandomAndSystematic', 'OnlySystematic'};

plotErrorNorms(names{1}, 10);
plotErrorNorms(names{1}, 5);
plotErrorNorms(names{1}, 0.5);

plotErrorNorms(names{2}, 10);
plotErrorNorms(names{2}, 5);
plotErrorNorms(names{2}, 0.5);

plotErrorNorms(names{3}, 10);
plotErrorNorms(names{3}, 5);
plotErrorNorms(names{3}, 0.5);
end

function plotErrorNorms(name, hmax)
load(['saveTmpErrorNormSimulations' name '.mat']);

kOrig= [0.0500
    0.1100
    0.2150
    1.2280
    0.2420
    0.0070];
lineStyle = '--';

for i = 1:len
    for j = 1:len
        p = norms(i);
        q = norms(j);
        
        subplot(len,len,(i-1)*len+j);
        [~,~,~,~,l] =plotKineticModelFit(models{i,j}.data.timeZ,models{i,j}.data.zml,models{i,j}.k,models{i,j}.z0opt);
        set(l,'Visible','off'); % no legend
        
        [t, z] = SolveKineticModel(models{i,j}.z0opt,kOrig);
        hold on
        h=plot(t,z(:,1),['b' lineStyle],t,z(:,2),['g' lineStyle],t,z(:,3),...
            ['r' lineStyle],t,z(:,4),['c' lineStyle],t,z(:,5),['m' lineStyle],...
            t,z(:,6),['k' lineStyle]);
        hold off
        axis([-Inf Inf 0 hmax])
        
        title(['L^p = ' num2str(p) ', L^q = ' num2str(q)]);
        drawnow;
    end
end

ViewPlot.Save(gcf,['Lpq-' name '-hmax' num2str(hmax)]);

end
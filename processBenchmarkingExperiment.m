function processBenchmarkingExperiment

good = processRunUniqunessExperimentRepetetiveFits();
% processUniqunessExperimentImportanceSampling()

end

function good = processRunUniqunessExperimentRepetetiveFits()
%100 reps as in the paper draft
load('Results\save_2013-08-10_125425BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_nonregularized_Noureddini');
load('Results\save_2013-08-06_171029BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_regularized_Noureddini');
%load('Results\save_2013-08-06_171301BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_regularized_Jansri');
load('Results\save_2013-08-12_085251BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_nonregularized_Jansri');
load('Results\save_2013-08-25_082739BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_onlyRandomError_nonregularized');
load('Results\save_2013-08-26_142908BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_onlyRandomError_nonregularized');

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
nvar
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
function processBenchmarkingExperiment

savefilenames = {'../data/save_2013-08-25_082739BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_onlyRandomError_nonregularized', ...
    '../data/save_2013-08-21_082412BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_onlyRandomError_regularized', ...
    '../data/save_2013-08-10_125425BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_nonregularized_Noureddini', ...
    '../data/save_2013-08-06_171029BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_regularized_Noureddini', ...
    '../data/save_2013-08-26_142908BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_onlyRandomError_nonregularized', ...
    '../data/save_2013-08-21_093754BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_onlyRandomError_regularized', ...
    '../data/save_2013-08-12_085251BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_nonregularized_Jansri', ...
    '../data/save_2013-08-06_171301BenchmarkingExperiment_100_Repetetive_Fits_RelativeLambda_regularized_Jansri', ...
    '../data/save_2020-05-12_232318BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_randomErr_nonregularized_minErr_nonregularized_minErr_001', ...
    '../data/save_2020-05-15_012451BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_randomErr_regularized_minErr_regularized_minErr_001', ...
    '../data/save_2020-05-16_135535BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_sysAndRandomErr_nonregularized_minErr_nonregularized_minErr_001', ...
    '../data/save_2020-05-18_163716BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Noureddini_sysAndRandomErr_regularized_minErr_regularized_minErr_001',...
    '../data/save_2020-05-25_015412BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_sysAndRandomErr_nonregularized_minErr_nonregularized_minErr_001.mat', ...
    '../data/save_2020-05-27_070850BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_sysAndRandomErr_regularized_minErr_regularized_minErr_001.mat', ... % _partial
    '../data/save_2020-05-23_023907BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_randomErr_regularized_minErr_regularized_minErr_001.mat', ...
    '../data/save_2020-05-20_160454BenchmarkingExperiment_100_Repetetive_Fits_NaN_RelativeLambda_Jansri_randomErr_nonregularized_minErr_nonregularized_minErr_001.mat'};

[good, allResults] = processLossRegResults(savefilenames);
[processed_all_results, tableString] = processAllResults(allResults);
saveTableAsTexFile(tableString);
end

function [allResults, tableString] = processAllResults(allResults)
alpha = 0.05;

allIds = allResults(:,1);
%  [allResults; problemId, problemNo, problemRand, problemSys, problemMin, ...
%                 regularized, i, j, kEuclideanMean(i,j), kEuclideanStd(i,j), sampleN]
isBest = zeros(size(allResults, 1), 1);
allResults = [allResults, isBest];
ids = unique(allResults(:,1));
for l = 1:(length(ids))
    id = ids(l);
    inds = allResults(:,1) == id;
    indsNum = 1:length(inds);
    indsNum = indsNum(inds);
    means = allResults(inds, 9);
    stds = allResults(inds, 10);
    ns = allResults(inds, 11);
    [minMean, indMin] = min(means);
    minSd = stds(indMin);
    minN = ns(indMin);
    for i = 1:sum(inds)
        alphaBonferroni = alpha / sum(inds);
        [hRejected, p] = my_ttest(minMean, means(i), minSd, stds(i), minN, ns(i), alphaBonferroni);
        if ~hRejected
            allResults(indsNum(i), 12) = 1;
        end
    end
end

    function pn = getReadablePnorms(pnorm)
        switch(pnorm)
            case 0.5
                pn = 'Root';
            case 1
                pn = 'Absolute';
            case 2
                pn = 'Square';
            case 'log'
                pn = 'Log';
            case 'rel'
                pn = 'Relative';
            otherwise
                error('processBenchmarkingExperiment:getReadablePnorms', 'Unsupported pnorm');
        end
    end

tableString = ['\\hline Exp. & Rand. & Sys. & Min. & Reg. & $L_\\text{time}$ & NaN & log & 0.5 & 1 & 2 & Wins '];


qnorms= {NaN, 'log', 0.5, 1, 2};
pnorms = {'log', 0.5, 1, 2};
plen = length(pnorms);
qlen = length(qnorms);
all_wins = zeros(plen, qlen);
for l = 1:(length(ids))
    id = ids(l);
    inds = allResults(:,1) == id;
    res = allResults(inds, :);
    regs = unique(allResults(inds, 6)).';
    if (res(1,2) == 1)
        experiment = 'N.';
    elseif (res(1,2) == 2)
        experiment = 'J.';
    elseif (res(1,2) == 3)
        experiment = 'K.';
    elseif (res(1,2) == 4)
        experiment = 'O.';
    else
        experiment = 'NA';
    end
    newId = ['\n ' '\\hline \\multirow{8}{*}{' experiment '} ' num2str(res(1,3:5), ' & \\\\multirow{8}{*}{$%d$} & \\\\multirow{8}{*}{$%d$} & \\\\multirow{8}{*}{$%d$}')];
    wins = zeros(1, qlen);
    for regularized = regs
        resr = res(res(:,6) == regularized, :);
        if (regularized == regs(1))
            newId = [newId '& \\multirow{4}{*}{$' num2str(regularized) '$}'];
        else
            newId = [' \\cline{5-12} \n & & & & \\multirow{4}{*}{$'  num2str(regularized)  '$}'];
        end
        
        for i = 1:plen
            resrp = resr(resr(:,7) == i, :);
            winsP = 0;
            tableString = [tableString ' \\\\ ' newId  ' & ' getReadablePnorms(pnorms{i}) '  '];
            newId = '\n & & & &';
            for j = 1:qlen
                resrpq = resrp(resrp(:,8) == j, :);
                kMean = resrpq(9);
                kStd = resrpq(10);
                mathbfOpen = '';
                mathbfClose = '';
                isWinner = resrpq(12) == 1;
                if (isWinner)
                    mathbfOpen = '\\boldsymbol{';
                    mathbfClose = '}';
                    wins(j) = wins(j) + 1;
                    winsP = winsP + 1;
                    all_wins(i,j) = all_wins(i,j) + 1;
                end
                
                tableString = [tableString ' & ' '\\cellcolor{green!' getCellColor(kMean) '}' ...
                    ' $' mathbfOpen num2str(kMean, '%2.1f') ' \\pm ' num2str(kStd, '%2.1f') mathbfClose '$ '];
                
            end
            tableString = [tableString ' & $' num2str(winsP) '$ '];
        end
    end
    tableString = [tableString ' \\\\ \\cline{5-12} \n \\multicolumn{6}{r}{$L_\\text{comp}$ wins} & ' num2str(wins, '$%d$ & $%d$ & $%d$ & $%d$ & $%d$') ' & '];
    
end
tableString = [tableString ' \\\\ \\hline \n'];
sprintf(tableString)
% Funciton processAllResults prints the same data as in tableString but
% sorted and with \boldfont added for the best results and those that are
% not statistically significantly worse than the best

% The number of wins for each method (for regularized regression only)
all_wins;
end


function [hRejected, p] = my_ttest(m1, m2, sd1, sd2, n1, n2, alpha)
if (nargin == 0)
    % Test case, expected output t = 1.959, dof = 7.031, (two-tailed) p =
    % 0.09077 % https://en.wikipedia.org/wiki/Student%27s_t-test#Unequal_variances
    A1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99];
    A2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98];
    m1 = mean(A1);
    m2 = mean(A2);
    sd1 = std(A1);
    sd2 = std(A2);
    n1 = length(A1);
    n2 = length(A2);
end
    function [x, y] = swap(x, y)
        z = x;
        x = y;
        y = z;
    end
if (m1 > m2)
    warning('my_ttest:Wrong order for paired test');
    [m1, m2] = swap(m1, m2);
    [sd1, sd2] = swap(sd1, sd2);
    [n1, n2] = swap(n1, n2);
end

s1n = (sd1^2)/n1;
s2n = (sd2^2)/n2;

% t statistic of the Welsch test (we have equal sample sizes)
t = (m1 - m2)/sqrt(s1n + s2n);

% degrees of freeedom

dof = ((s1n + s2n)^2) / ((s1n^2)/(n1-1) + (s2n^2)/(n2-1));

p = tcdf(t, dof); % one-sided p-value
hRejected = p <= alpha;
end

function [text, id, idNumeric, regularized, ...
    problemNo, problemRand, problemSys, problemMin] = extract_design_from_filename(filename)
load(filename)

if (all(sysErr == 0))
    sysErr = false;
else
    sysErr = true;
end
minErrStr = false;
if (exist('minErr', 'var'))
    if (minErr > 0)
        minErrStr = true;
    end
end

regularized = NaN;
if (~isempty(regexp(filename, '_nonregularized', 'once')))
    regularized = false;
end
if (~isempty(regexp(filename, '_regularized', 'once')))
    regularized = true;
end

experiment = NaN;

if strcmp(id,'Noureddini')
    experiment = 'N.';
    experimentNumeric = 1;
end
if strcmp(id,'Jansri')
    experiment = 'J.';
    experimentNumeric = 2;
end
if strcmp(id,'Klofutar')
    experiment = 'K.';
    experimentNumeric = 3;
end
if strcmp(id,'Oh')
    experiment = 'O.';
    experimentNumeric = 4;
end

text = [experiment ' & ' '$1$' ' & $' num2str(sysErr) '$ & $' num2str(minErrStr) '$ & $'  num2str(regularized) '$ '];
id = [experiment '_1_' num2str(sysErr) '_' num2str(minErrStr)];
idNumeric = experimentNumeric * 1e6 + 1e5 + sysErr*1e4 + minErrStr * 1e3;
problemNo = experimentNumeric;
problemRand = 1;
problemSys = sysErr * 1;
problemMin = minErrStr * 1;
end

function col = getCellColor(k)
mmax = 50;
revk = max(0,(mmax - k))* (100/mmax);
col = round(100 * (revk / 100).^1 );
col(col <= 95) = 50/95 * col(col<=95);
col(col>95) = -900 + 10 * col(col>95);
%plot(k, col)
col = num2str(round(col));
end

function [good, allResults] = processLossRegResults(savefilenames)

tableString = '';
plotting = true;
allResults = [];

allP = [];

for indFile = 1:length(savefilenames)
    load(savefilenames{indFile});
    fprintf('--------------------------------\n')
    disp(['File: ' savefilenames{indFile} '\n']);
    if(plotting)
        figure()
    end
    
    if (isempty(tableString))
        tableString = ['\\hline Exp. & Rand. & Sys. & Min. & Reg. & $L_\\text{time}$ & ' ...
            num2str(qnorms{1}) ' & ' num2str(qnorms{2}) ' & ' num2str(qnorms{3}) ' & ' num2str(qnorms{4}) ' & ' num2str(qnorms{5}) ' '];
    end
    
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
                if isempty(models{r,i,j})
                    mtr(r,:) = NaN(1,6);
                else
                    mtr(r,:) = models{r,i,j}.k.';
                end
            end
            kVal{i,j}= mtr;
        end
    end
    good=NaN;
    kEuclideanMean = NaN;
    kEuclideanStd = NaN;
    testId = 1;
    dunnErrors = NaN(plen*qlen*rep,2);
    bestId = [1, 1, 1];
    bestIdOneStep = [1,1];
    [designText, designId, problemId, regularized, ...
        problemNo, problemRand, problemSys, problemMin] = extract_design_from_filename(savefilenames{indFile});
    for i = 1:plen
        tableString = [tableString ' \\\\ \n  ' designText ' & ' num2str(pnorms{i}) '  '];
        for j = 1:qlen
            kf = sum(kVal{i,j}(:,[1 3 5]).').';
            kb = sum(kVal{i,j}(:,[2 4 6]).').';
            
            
            
            
            [r, ~] = size(kVal{i,j});
            isGood = all(kVal{i,j}>repmat(data.k.'*lRelativeAccuracy,r,1),2) & ...
                all(kVal{i,j}<repmat(data.k.'*uRelativeAccuracy,r,1),2);
            kc = kVal{i,j}(isGood,:);
            
            good(i,j) = round(100*numel(kc)/numel(kVal{i,j}));
            
            if(plotting)
                subplot(plen,qlen,qlen*(i-1)+j)
                h = plot(kf,kb,'go',sum(data.k([1,3,5])),sum(data.k([2,4,6])),'*');
                set(gca,'xScale','log','yScale','log');
                axis([1e-2 1e3 1e-2 1e3]);
                xlabel('k_f = k_1 + k_3 +k_5')
                ylabel('k_b = k_2 + k_4 +k_6')
                hl = line(sum(kc(:,[1 3 5]).').',sum(kc(:,[2 4 6]).').');
                set(hl,'color','r','marker','o','linestyle','none');
            end
            
            kActualMatrix = repmat(data.k.', size(kVal{i,j}, 1),1);
            kEuclideanDists = sqrt(sum(((kVal{i,j} - kActualMatrix).^2).')).';
            kEuclideanMean(i,j) = mean(kEuclideanDists);
            kEuclideanStd(i,j) = std(kEuclideanDists);
            nvar(i,j) = norm(kVal{i,j});
            if all(~isnan(kEuclideanDists))
                [~,P,~] = swtest(kEuclideanDists);
                allP = [allP; P kEuclideanMean(i,j)];
            end
            
            sampleN = length(kEuclideanDists);
            allResults = [allResults; problemId, problemNo, problemRand, problemSys, problemMin, ...
                regularized, i, j, kEuclideanMean(i,j), kEuclideanStd(i,j), sampleN];
            
            tableString = [tableString ' & ' '\\cellcolor{green!'	getCellColor(kEuclideanMean(i,j))	'}' ...
                ' $' num2str(kEuclideanMean(i,j), '%2.1f') ' \\pm ' num2str(kEuclideanStd(i,j), '%2.1f') '$ '];
            
            if(plotting)
                h = boxplot(kVal{i,j});
                axis([0 7 1e-3 1e3]);
                set(gca,'yScale','log');
                hl = line([1:6], data.k, 'MarkerSize', 10, 'LineStyle', ':', 'Marker', 'o', 'LineWidth', 2);
                
                title(['p = ' num2str(pnorms{i}) ', q = ' num2str(qnorms{j})]);
                xlabel('Index i')
                ylabel('Rate constant k_i')
                title(['p = ' num2str(pnorms{i}) ', q = ' num2str(qnorms{j}),... % ', good = ', num2str(good(i,j)) '%']);
                    ', ||k*-k|| = ', sprintf('%1.1f', kEuclideanMean(i,j)), ' � ', sprintf('%1.1f', kEuclideanStd(i,j))]);
            end
            if(j == 1)
                if(good(i,1) >= good(1,bestIdOneStep(1)))
                    bestIdOneStep(1) = i;
                    bestIdOneStep(2) = testId;
                end
            end
            if (good(i,j) >= good(bestId(1),bestId(2)))
                bestId(1) = i;
                bestId(2) = j;
                bestId(3) = testId;
            end
            kErr=kVal{i,j} - repmat(data.k.',r,1); % errors of each rate constant
            kErrNorm = sum((kErr.').^2).';
            dunnErrors((1+(testId-1)*rep):(testId*rep),:) = ...
                [isGood, repmat(testId,rep,1)];
            %             [kErrNorm, repmat(testId,rep,1)];
            testId = testId+1;
        end
    end
    % good
    %nvar
    %round(nvar)
    %figure()
    %bar3(nvar)
    %bar3(good)
    control = dunnErrors(dunnErrors(:,2) == bestId(3),1);
    for l = 1:plen*qlen
        group = dunnErrors(dunnErrors(:,2) == l,1);
        p(l) = signtest(control,group,'tail','right');
    end
    h = p<0.05/plen*(qlen-1);
    
    controlOneStep = dunnErrors(dunnErrors(:,2) == bestIdOneStep(2),1);
    for l = 1:plen
        ind = [1 6 11 16];
        groupOneStep = dunnErrors(dunnErrors(:,2) == ind(l),1);
        pOneStep(l) = signtest(controlOneStep,groupOneStep,'tail','right');
    end
    hOneStep = pOneStep<0.05/8;
    
    
    testId = 1;
    pPrint = '';
    hPrint = '';
    for i = 1:plen
        for j = 1:qlen
            pPrint = [pPrint sprintf([' ' num2str(p(testId)) ' '])];
            if (h(testId))
                htmp = 1;
            else
                htmp = 0;
            end
            hPrint = [hPrint sprintf([' ' num2str(htmp) ' '])];
            testId = testId+1;
            
            if (j == 1)
                pPrint = [pPrint '   '];
                hPrint = [hPrint '   '];
            end
        end
        
        pPrint = [pPrint '\n'];
        hPrint = [hPrint '\n'];
    end
    % fprintf(pPrint)
    
    %     fprintf(hPrint)
    %     kEuclideanMean
    %     good
    %     hOneStep
    %     pOneStep
end
printingFinalResults = false;
if (printingFinalResults)
    tableString = [tableString '\\\\ \\hline \n'];
    fprintf(tableString)
    allP(allP(:,1)>0.05,:) % Results of Shapiro-Wilk normality tests
end
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

function saveTableAsTexFile(tableString)
preText = {'\\documentclass{article}',...
'\\usepackage{amssymb}',...
'\\usepackage{amsmath}',...
'\\usepackage{rotating}',...
'\\usepackage{multirow}',...
'\\usepackage[table]{xcolor}',...
'\\usepackage{tabularx}',...
'\\begin{document}',...
'\\begin{table}[htbp]',...
'\\caption{Mean $\\pm$ standard deviation of Euclidean distance between the optimized $\\boldsymbol k$ and underlying $\\boldsymbol k^*$ rate constants for $100$ independent simulations \\label{tab:resultsAll}}',...
'\\tiny',...
'\\begin{tabular}{cccccccccccc}'};
postText = {'\\end{tabular}',...
'\\end{table}',...
'\\end{document}'};
allText = [sprintf('%s\n',preText{:}) tableString sprintf('%s\n', postText{:})];
fn = '../results/able_5_loss_reg_comparison.tex';
fid = fopen(fn,'w');
fprintf(fid,allText);
fclose(fid);
disp(['Saved table as ' fn]);
end
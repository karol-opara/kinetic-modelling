function processOptimizationComparison

[good, allResults] = processFiles();
% processUniqunessExperimentImportanceSampling()

% processed_all_results = processAllResults(allResults);

end

function allResults = processAllResults(allResults)
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
all_wins
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


function col = getCellColor(k)
mmax = 50;
revk = max(0,(mmax - k))* (100/mmax);
col = round(100 * (revk / 100).^1 );
col(col <= 95) = 50/95 * col(col<=95);
col(col>95) = -900 + 10 * col(col>95);
%plot(k, col)
col = num2str(round(col));
end

function [good, allResults] = processFiles()
savefilenames = {'Results/save_2021-08-18_201142_OptimizationComparison_30_RepetetiveFits_4-optimizers_nonregularized_minErr_001.mat',...
    'Results/save_2021-08-19_043100_OptimizationComparison_30_RepetetiveFits_4-optimizers_nonregularized_minErr_001.mat'};% {'Results\save_2021-08-15_112810OptimizationComparison_100_RepetetiveFits_test_nonregularized_minErr_001.mat'};

if false == exist('pqnames')
    pqnames = {'Relative', 'Square', 'Regularized log-square'};
end

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
        tableString = ['\\hline Exp. & $L_\\text{time}$' ];
        for io = 1:olen
            tableString = [tableString ' & ' optimizers{io}];
        end
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
    
    
    kVal = cell(plen,olen);
    for j = 1:plen
        for i = 1:olen
            mtr = NaN(rep,6);
            for (r = 1:rep)
                if isempty(models{r,j,i})
                    mtr(r,:) = NaN(1,6);
                else
                    mtr(r,:) = models{r,j,i}.k.';
                end
            end
            kVal{j,i}= mtr;
        end
    end
    good=NaN;
    kEuclideanMean = NaN;
    kEuclideanStd = NaN;
    testId = 1;
    dunnErrors = NaN(plen*olen*rep,2);
    bestId = [1, 1, 1];
    bestIdOneStep = [1,1];
    for i = 1:plen
        tableString = [tableString ' \\\\ \n  ' id ' & ' pqnames{i} '  '];
        for j = 1:olen
            kf = sum(kVal{i,j}(:,[1 3 5]).').';
            kb = sum(kVal{i,j}(:,[2 4 6]).').';
            

            [r, ~] = size(kVal{i,j});
            isGood = all(kVal{i,j}>repmat(data.k.'*lRelativeAccuracy,r,1),2) & ...
                all(kVal{i,j}<repmat(data.k.'*uRelativeAccuracy,r,1),2);
            kc = kVal{i,j}(isGood,:);
            
            good(i,j) = round(100*numel(kc)/numel(kVal{i,j}));
            
            if(plotting)
                disp(['Number ', num2str(olen*(i-1)+j)])
                        subplot(plen,olen,olen*(i-1)+j)
                        h = plot(kf,kb,'go',sum(data.k([1,3,5])),sum(data.k([2,4,6])),'*');
                        set(gca,'xScale','log','yScale','log');
                        axis([1e-2 1e3 1e-2 1e3]);
                        xlabel('k_f = k_1 + k_3 +k_5')
                        ylabel('k_b = k_2 + k_4 +k_6')
                        hl = line(sum(kc(:,[1 3 5]).').',sum(kc(:,[2 4 6]).').');
                        set(hl,'color','r','marker','o','linestyle','none');
                        title([pqnames{i} ' ' optimizers{j}]);
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
%             allResults = [allResults; id, pqnames{i}, optimizers{j}, kEuclideanMean(i,j), kEuclideanStd(i,j), sampleN];
            
            tableString = [tableString ' & ' '\\cellcolor{green!'	getCellColor(kEuclideanMean(i,j))	'}' ...
                ' $' num2str(kEuclideanMean(i,j), '%2.1f') ' \\pm ' num2str(kEuclideanStd(i,j), '%2.1f') '$ '];
            
            if(plotting)
                h = boxplot(kVal{i,j});
                axis([0 7 1e-3 1e3]);
                set(gca,'yScale','log');
                hl = line([1:6], data.k, 'MarkerSize', 10, 'LineStyle', ':', 'Marker', 'o', 'LineWidth', 2);
                
                xlabel('Index i')
                ylabel('Rate constant k_i')
                title([pqnames{i} ' ' optimizers{j} ...
                    ', ||k*-k|| = ', sprintf('%1.1f', kEuclideanMean(i,j)), ' ± ', sprintf('%1.1f', kEuclideanStd(i,j))]);
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
    nvar
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
    
    fprintf(hPrint)
    kEuclideanMean
    good
    hOneStep
    pOneStep
end
tableString = [tableString '\\\\ \\hline \n'];
fprintf(tableString)
allP(allP(:,1)>0.05,:) % Results of Shapiro-Wilk normality tests 
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
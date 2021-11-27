function processOptimizationComparison
    
savefilenames = {'../data/regularization/save_2021-09-07_161504_OptimizationComparison_Noureddini_30_RepetetiveFits_7-optimizers_nonregularized_minErr_001_pnameFixed.mat',...
    '../data/regularization/save_2021-09-08_024025_OptimizationComparison_Jansri_30_RepetetiveFits_7-optimizers_nonregularized_minErr_001_pnameFixed.mat',...
    '../data/regularization/save_2021-09-08_154749_OptimizationComparison_Klofutar_30_RepetetiveFits_7-optimizers_nonregularized_minErr_001_pnameFixed.mat',...
    '../data/regularization/save_2021-09-09_121801_OptimizationComparison_Noureddini_30_RepetetiveFits_7-optimizers_regularized_minErr_001.mat',...
    '../data/regularization/save_2021-09-09_151122_OptimizationComparison_Jansri_30_RepetetiveFits_7-optimizers_regularized_minErr_001.mat',...
    '../data/regularization/save_2021-09-09_191054_OptimizationComparison_Klofutar_30_RepetetiveFits_7-optimizers_regularized_minErr_001.mat'};

[good, allResults] = processFiles(savefilenames);
[processed_all_results, tableString]  = processAllResults(allResults);
saveTableAsTexFile(tableString);
end

function [allResults, tableString] = processAllResults(allResults)
alpha = 0.05;

% allResults{ar_row,1} = id;
% allResults{ar_row,2} = pnames{i};
% allResults{ar_row,3} = optimizers{j};
% allResults{ar_row,4} = kEuclideanMean(i,j);
% allResults{ar_row,5} = kEuclideanStd(i,j);
% allResults{ar_row,6} = sampleN;
% allResults{ar_row,7} = 0; % is_best -- to be filled in later
% allResults{ar_row,8} = [id, '-', optimizers{j}]; % problem-algorithm pair

ids = unique(allResults(:,8));
for l = 1:(length(ids))
    id = ids{l};
    inds = strcmp(allResults(:,8), id);
    indsNum = 1:length(inds);
    indsNum = indsNum(inds);
    means = cell2mat(allResults(inds, 4));
    stds = cell2mat(allResults(inds, 5));
    ns = cell2mat(allResults(inds, 6));
    [minMean, indMin] = min(means);
    minSd = stds(indMin);
    minN = ns(indMin);
    for i = 1:sum(inds)
        alphaBonferroni = alpha / sum(inds);
        [hRejected, p] = my_ttest(minMean, means(i), minSd, stds(i), minN, ns(i), alphaBonferroni);
        if ~hRejected
            allResults(indsNum(i), 7) = {1};
        end
    end
end


ids = unique(allResults(:,1));
problems = unique(allResults(:,2));
plen = length(problems);
optimizers = unique(allResults(:,3));
olen = length(optimizers);
all_wins = zeros(plen, olen);

tableString = ['\\hline Exp. & Formulation'];
for i = 1:olen
    tableString = [tableString ' & ', makeReadableOptimizer(optimizers{i})];
end
tableString = [tableString  ' & Wins '];

wins = zeros(1,olen);
for l = 1:(length(ids))
    id = ids{l};
    inds = strcmp(allResults(:,1), id);
    indsNum = 1:length(inds);
    indsNum = indsNum(inds);
    res = allResults(indsNum, :);
    
    if (strcmp(id, 'Noureddini'))
        experiment = 'N.';
    elseif (strcmp(id, 'Jansri'))
        experiment = 'J.';
    elseif (strcmp(id, 'Klofutar'))
        experiment = 'K.';
    else
        experiment = 'NA';
    end
    
    hline = true;
    
    for i = 1:plen
        resp = res(strcmp(res(:,2), problems{i}), :);
        winsP = 0;
        
        pname = problems{i};
        if (strcmp(problems{i}, 'Log-square'))
            pname = 'Log-sq.';
        end
        if (strcmp(problems{i}, 'Regularized log-square'))
            pname = 'Reg. log-sq.';
        end
        
        hl = '';
        if hline
            hl = ['\\hline \n \\multirow{4}{*}{' experiment '} '];
            hline = false;
        end
        
        tableString = [tableString ' \\\\ \n  ' hl ' & ' pname '  '];
        
        for j = 1:olen
            respq = resp(strcmp(resp(:,3), optimizers{j}), :);
            kMean = cell2mat(respq(4));
            kStd = cell2mat(respq(5));
            mathbfOpen = '';
            mathbfClose = '';
            isWinner = cell2mat(respq(7)) == 1;
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
compareWins = false;
if (compareWins)
    tableString = [tableString ' \\\\ \\hline \n & Wins '];
    for j = 1:olen
        tableString = [tableString ' & ', num2str(wins(j))];
    end
    tableString = [tableString ' & \\\\ \\hline \n'];
else
    tableString = [tableString ' \\\\ \\hline \n'];
end

sprintf(tableString)
% Funciton processAllResults prints the same data as in tableString but
% sorted and with \boldfont added for the best results and those that are
% not statistically significantly worse than the best

% The number of wins for each method (for regularized regression only)
all_wins
sum(all_wins.').'
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

function opt = makeReadableOptimizer(str)
opt = str;
switch str %make it case insensitive. all lower cases
    case 'ampso'
        opt = 'AMPSO';
    case 'apgskimode'
        opt = 'APGSK-IMODE';
    case 'cmaes'
        opt = 'CMA-ES';
    case 'derandinfty'
        opt = 'DE/rand/$\\infty$';
    case 'fmincon'
        opt = 'SQP';
    case 'madDE'
        opt = 'MadDE';
    case 'somat3a'
        opt = 'SOMA T3A';
    otherwise
        error('processOptimizationComparison:Not supported optimizer type');
end
end

function [good, allResults] = processFiles(savefilenames)

% if false == exist('pqnames')
%     pqnames = {'Relative', 'Square', 'Regularized log-square'};
% end

tableString = '';
plotting = true;
allResults = cell(1,7); % plen * olen * experiments
ar_row = 0;

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
        tableString = [tableString ' \\\\ \n  ' id ' & ' pnames{i} '  '];
        for j = 1:olen
            kf = sum(kVal{i,j}(:,[1 3 5]).').';
            kb = sum(kVal{i,j}(:,[2 4 6]).').';
            
            data.k = dataN(1).k;
            
            [r, ~] = size(kVal{i,j});
            isGood = all(kVal{i,j}>repmat(data.k.'*lRelativeAccuracy,r,1),2) & ...
                all(kVal{i,j}<repmat(data.k.'*uRelativeAccuracy,r,1),2);
            kc = kVal{i,j}(isGood,:);
            
            good(i,j) = round(100*numel(kc)/numel(kVal{i,j}));
            
            if(plotting)
                % disp(['Number ', num2str(olen*(i-1)+j)])
                subplot(plen,olen,olen*(i-1)+j)
                h = plot(kf,kb,'go',sum(data.k([1,3,5])),sum(data.k([2,4,6])),'*');
                set(gca,'xScale','log','yScale','log');
                axis([1e-2 1e3 1e-2 1e3]);
                xlabel('k_f = k_1 + k_3 +k_5')
                ylabel('k_b = k_2 + k_4 +k_6')
                hl = line(sum(kc(:,[1 3 5]).').',sum(kc(:,[2 4 6]).').');
                set(hl,'color','r','marker','o','linestyle','none');
                title([pnames{i} ' ' optimizers{j}]);
            end
            
            kActualMatrix = repmat(data.k.', size(kVal{i,j}, 1),1);
            kEuclideanDists = sqrt(sum(((kVal{i,j} - kActualMatrix).^2).')).';
            kEuclideanMean(i,j) = mean(kEuclideanDists);
            kEuclideanStd(i,j) = std(kEuclideanDists);
            nvar(i,j) = norm(kVal{i,j});
            if all(~isnan(kEuclideanDists))
                if length(kEuclideanDists) < 3
                    P = NaN;
                else
                    [~,P,~] = swtest(kEuclideanDists);
                end
                allP = [allP; P kEuclideanMean(i,j)];
            end
            
            sampleN = length(kEuclideanDists);
            ar_row = ar_row + 1;
            allResults{ar_row,1} = id;
            allResults{ar_row,2} = pnames{i};
            allResults{ar_row,3} = optimizers{j};
            allResults{ar_row,4} = kEuclideanMean(i,j);
            allResults{ar_row,5} = kEuclideanStd(i,j);
            allResults{ar_row,6} = sampleN;
            allResults{ar_row,7} = 0; % is_best -- to be filled in later
            allResults{ar_row,8} = [id, '-', optimizers{j}]; % problem-algorithm pair
            
            tableString = [tableString ' & ' '\\cellcolor{green!'	getCellColor(kEuclideanMean(i,j))	'}' ...
                ' $' num2str(kEuclideanMean(i,j), '%2.1f') ' \\pm ' num2str(kEuclideanStd(i,j), '%2.1f') '$ '];
            
            if(plotting)
                h = boxplot(kVal{i,j});
                axis([0 7 1e-3 1e3]);
                set(gca,'yScale','log');
                hl = line([1:6], data.k, 'MarkerSize', 10, 'LineStyle', ':', 'Marker', 'o', 'LineWidth', 2);
                
                xlabel('Index i')
                ylabel('Rate constant k_i')
                title([pnames{i} ' ' optimizers{j} ...
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
    for l = 1:plen*olen
        group = dunnErrors(dunnErrors(:,2) == l,1);
        p(l) = signtest(control,group,'tail','right');
    end
    h = p<0.05/plen*(olen-1);
    
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
        for j = 1:olen
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
    
    printInterimResults = false;
    if (printInterimResults)
        fprintf(hPrint)
        kEuclideanMean
        good
        hOneStep
        pOneStep
    end
end
printFinalResults = false;
if (printFinalResults)
tableString = [tableString '\\\\ \\hline \n'];
fprintf(tableString)
allP(allP(:,1)>0.05,:) % Results of Shapiro-Wilk normality tests
end
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
'\\begin{table}[h]',...
'    \\centering',...
'    \\tiny',...
'    \\caption{Errors obtained for three problem formulations by seven optimization algorithms for simulating three experiments}',...
'    \\label{tab:comparisonProblemsAlgorithms}',...
'    \\begin{tabular}{cccccccccc}'};
postText = {'\\end{tabular}',...
'\\end{table}',...
'\\end{document}'};
allText = [sprintf('%s\n',preText{:}) tableString sprintf('%s\n', postText{:})];
fn = '../results/Table_7_optimizers_comparison.tex';
fid = fopen(fn,'w');
fprintf(fid,allText);
fclose(fid);
disp(['Saved table as ' fn]);
end
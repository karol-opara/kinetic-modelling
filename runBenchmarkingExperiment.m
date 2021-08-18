function runBenchmarkingExperiment(name, rndErr, sysErr, minErr, dataset, lambda, N, lambdas, compareOptimizers)

if (nargin < 2)
    rndErr = 0.15;
    sysErr = [1 0 0 0 0 0];
end
if (nargin < 5)
    dataset = 'Oh';
    warning('runBenchmarkingExperiment:NoExperimentID','No experiment ID provided');
end

% if (areTheseToolboxesInstalled({'MATLAB','Parallel Computing Toolbox'}))
%     matlabpool(feature('numCores')-1);
% end

% RunUniqunessExperimentImportanceSampling(randErr, sysErr)
if compareOptimizers
    RunOptimizerComparison(name, rndErr, sysErr, minErr, dataset, N, lambdas);
else
    RunUniqunessExperimentRepetetiveFits(name, rndErr, sysErr, minErr, dataset,lambda,N,lambdas);
end


end


function RunUniqunessExperimentRepetetiveFits(name, rndErr, sysErr, minErr, id,lambda, N,lambdas)
savefilename = ['Results/' 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    'BenchmarkingExperiment_' num2str(N) '_Repetetive_Fits_' num2str(lambda)...
    '_RelativeLambda_' name '_minErr_' num2str(minErr)];
savefilename(ismember(savefilename,' ,.:;!'))=[];

%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only 10 repeats');
pnorms = {'log', 0.5, 1, 2};
qnorms = {NaN,'log', 0.5, 1, 2};
%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only NaN norms tried');
plen = length(pnorms);
qlen = length(qnorms);
models = cell(N,plen,plen);
for i = 1:N
    dataN(i) = CreateBenchmarkProblem(rndErr, sysErr, minErr, id);
end
data = dataN(1);
for i = 1:plen
    %fprintf(['\n' num2str(rep) ': ']);
    for j = 1:qlen
        lambda = [1 lambdas(i,j)];
        parfor rep = 1:N
            data = dataN(rep);
            p = pnorms{i};
            q = qnorms{j};
            
            models{rep,i,j} = EstimateKineticModel(data,p,q,'batch',lambda*lossFunctionOrderMultiplier(i,j),...
                struct(), 'cmaes');
            
            %         subplot(len,len,(i-1)*len+j);
            %         plotKineticModelFit(models{i,j}.data.timeZ,models{i,j}.data.zml,models{i,j}.k,models{i,j}.z0opt);
            %         title(['L^p = ' num2str(p) ', L^q = ' num2str(q)]);
            %         drawnow;
        end
        fprintf('.');
        save([savefilename '_partial']);
    end
end
save(savefilename);
disp(['Saved as ' savefilename])
end

function lossMult = lossFunctionOrderMultiplier(i,j)
load('Results/save_RegularizationCoefficients_2013-07-23_092626_1e5Dim_LBFGS_PoorData_NonregularizedErrors');
lossMult = errMin(i,j);
end

% function RunUniqunessExperimentImportanceSampling(rndErr, sysErr, name)
% savefilename = ['Results/' 'save_BenchmarkingExperiment_ImportanceSamping_' ...
%     datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];
% norms = [0.5];
% len = length(norms);
% models = cell(len,len);
% samples = cell(len,len);
%
% N = 1e2;
%
% data = CreateBenchmarkProblem(rndErr, sysErr, id);
% for i = 1:len
%     for j = 1:len
%         p = norms(i);
%         q = norms(j);
%
%         samples{i,j} = UniqunessOfKineticSolution(p,q,data,N);
%         %models{i,j} = EstimateKineticModel(data,p,q);
%         save(['saveTmpBenchmarkingExperiment_' name]);
%
%         %         subplot(len,len,(i-1)*len+j);
%         %         plotKineticModelFit(models{i,j}.data.timeZ,models{i,j}.data.zml,models{i,j}.k,models{i,j}.z0opt);
%         %         title(['L^p = ' num2str(p) ', L^q = ' num2str(q)]);
%         %         drawnow;
%     end
% end
% save(savefilename);
% end

% function [sample, iter] = UniqunessOfKineticSolution(p,q,data, N)
% % Use the importance sampling to assess uniquness of solution
% kExperimental = [
%     5.000E-02	1.100E-01	2.150E-01	1.228E+00	2.420E-01	7.000E-03
%     1.030E-01	3.100E-02	6.300E-02	1.000E-02	1.600E-02	1.750E-01
%     2.860E-02	1.440E-02	5.800E-03	2.130E-02	1.110E-02	5.000E-04
%     8.000E-01	5.950E+00	1.050E+01	1.590E+01	3.400E-01	3.500E-03
%     1.550E+00	8.500E+00	2.050E+01	2.250E+01	6.100E-01	1.200E-03
%     2.050E+00	1.090E+01	3.010E+01	2.950E+01	8.300E-01	1.000E-04
%     1.500E+00	1.370E+01	2.300E+01	4.140E+01	4.000E-01	2.600E-03
%     3.060E+00	2.390E+01	3.250E+01	5.750E+01	5.400E-01	9.000E-04
%     4.000E+00	2.700E+01	5.500E+01	6.550E+01	9.100E-01	1.000E-04
%     2.579E+00	2.000E-02	6.000E-01	1.010E-01	9.000E-01	2.100E-02
%     2.600E+00	2.480E-01	1.186E+00	2.270E-01	2.303E+00	2.200E-02
%     2.620E+00	7.000E-01	1.210E+00	4.000E-01	2.360E+00	2.800E-02
%     2.470E-02	5.580E-02	7.030E-02	1.500E-03	3.940E-02	4.200E-03
%     7.720E-02	1.680E-01	9.720E-02	2.650E-02	6.700E-02	8.800E-03
%     4.430E-02	2.334E-01	6.450E-02	6.990E-02	2.681E-01	4.700E-03
%     8.790E-02	4.777E-01	1.555E-01	1.396E-01	7.478E-01	6.100E-03
%     1.057E-02	0.000E+00	1.184E-01	8.187E-02	1.310E-01	2.010E-03
%     ];
%
% kMean = 1*mean(kExperimental);
%
% k = NaN(2*N,6);
% logW = NaN(2*N,1);
% i = 0;
% iter = 0;
%
% for i = 1:(10*N)
%     kRnd=exprnd(rand()*kMean).';
%     logQ(i) = sum(log(kRnd)) - (1./kMean) * kRnd; % logarithm of proposal density
%     err = ObjectiveFunction(kRnd, data.zml, data.timeZ, data.z0ml,p,q);
%     logAP(i) = log(1/err); % logarithm of an unnormalized true density
%
%     logW(i) = (logAP(i) - logQ(i));
%     k(i,:) = kRnd;
% end
%
% w = exp(logW);
% w = w/sum(w);
%
% ind = randp(w,N,1);
% sample = k(ind,:);
%
% figure()
% set(gcf,'Position',[100 100 1000 450]);
% subplot(1,2,1)
% gini = 1 - sum(w.^2);
% plot(cumsum(w));
% title(['L^p = ' num2str(p)...
%     ', L^q = ' num2str(q) ', Gini = ' sprintf('%0.2f',gini)]);
% axis([-Inf Inf 0 1]);
% xlabel('Observation index');
% ylabel('CDF of weights');
%
% subplot(1,2,2)
% epsilon = 1e-1*(rand(N,2)-0.5);
% plot(sum(k(:,[1 3 5]).').',sum(k(:,[2 4 6]).').','g.',...
%     sum(sample(:,[1 3 5]).').'.*(1+epsilon(:,1)),...
%     sum(sample(:,[2 4 6]).').'.*(1+epsilon(:,2)),'bo');
% set(gca,'XScale','log','YScale','log');
% title([num2str(length(unique(ind))) ' unique out of '   ...
%     num2str(N) ' (' sprintf('%2.0f',100*length(unique(ind))/N) '%)']);
% axis([1e-1 1e2 1e-1 2e2]);
% xlabel('k_f = k_1 + k_3 + k_5');
% ylabel('k_b = k_2 + k_4 + k_6');
% legend('proposal','resampling');
% end

function data = CreateBenchmarkProblem(rndErr, sysErr, minErr, id)
% first
if (strcmp(id,'Oh'))
    k = [0.0500
        0.1100
        0.2150
        1.2280
        0.2420
        0.0070];
elseif (strcmp(id,'Jansri'))
    k = [2.6
        0.248
        1.186
        0.227
        2.303
        0.022];
elseif (strcmp(id,'Noureddini'))
    k = [0.0500
        0.1100
        0.2150
        1.2280
        0.2420
        0.0070];
elseif (strcmp(id,'Klofutar'))
    k = [0.0443
        0.2334
        0.0645
        0.0699
        0.2681
        0.0047];
else
    error('runBenchmarkingExperiment:WrongID','Wrong ID of the reaction rates');
end

conc = [6 0 1 0 0 0];

% totalMolePerLitre = 10; % in simulations and initial submission
totalMolePerLitre = 0.789359594 + 0.057795443 + 0.003075868 + 0.003075868 + 5.002450688; % in revision

z0ml = totalMolePerLitre * conc/sum(conc);
if (strcmp(id,'Oh'))
    timeZ = [2     4     6     8    10    12    20    30    40    50    60];
elseif(strcmp(id,'Jansri'))
    timeZ = [0.5, 1, 3, 5, 7, 9, 12, 15, 18, 20 ];
elseif (strcmp(id,'Noureddini'))
    timeZ = [1 2 3 4 5 6 8 10 15 20 30 50 60 90]; % Nouredini and Zhou (1997)
elseif (strcmp(id, 'Klofutar'))
    timeZ = [2 8 12 16 20 30 50 70 90]; % Pixel-wise reverse engineering from Fig. 4B
end
[t, z] = SolveKineticModel(z0ml,k);
zKinetic = interp1(t,z,timeZ);

%rng('shuffle');

rnd = rndErr* sqrt(pi/2) * randn(length(timeZ),length(z0ml)) + minErr * randn(length(timeZ),length(z0ml));
zml = zKinetic+zKinetic.*rnd+repmat(sysErr,length(timeZ),1);

data = struct('k',k,'z0ml',z0ml,'timeZ',timeZ,'zml',zml,'zKinetic',zKinetic);

end

function RunOptimizerComparison(name, rndErr, sysErr, minErr, id, N,lambdas)
savefilename = ['Results/' 'save_' datestr(now,'yyyy-mm-dd_HHMMSS') ...
    '_OptimizationComparison_' num2str(N) '_RepetetiveFits_' ...
    name '_minErr_' num2str(minErr)];
savefilename(ismember(savefilename,' ,.:;!'))=[];

%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only 10 repeats');
pnorms = {'rel', 2, 2};
qnorms = {NaN, NaN, 'log'};
pqnames = {'Relative', 'Square', 'Regularized log-square'};
optimizers = {'madDE', 'fmincon', 'cmaes', 'derandinfty'};
%warning('runBenchmarkingExperiment:RunUniqunessExperimentRepetetiveFits','Only NaN norms tried');
plen = length(pnorms);
qlen = length(qnorms);
olen = length(optimizers);
models = cell(N,olen,plen);
for i = 1:N
    dataN(i) = CreateBenchmarkProblem(rndErr, sysErr, minErr, id);
end
data = dataN(1);
for i = 1:plen
    %fprintf(['\n' num2str(rep) ': ']);
    lambda = [1 0];
    if(strcmp(qnorms{i},'log'))
        lambda = [1 lambdas(4,2)] * lossFunctionOrderMultiplier(4,2);
    end
    for j = 1:olen
        optimizer = optimizers{j};
        parfor rep = 1:N
        % for rep = 1:N
            data = dataN(rep);
            p = pnorms{i};
            q = qnorms{i};
            options = struct();
            
            models{rep,i,j} = EstimateKineticModel(data,p,q,'batch',...
                lambda, options, optimizer);
            
            %         subplot(len,len,(i-1)*len+j);
            %         plotKineticModelFit(models{i,j}.data.timeZ,models{i,j}.data.zml,models{i,j}.k,models{i,j}.z0opt);
            %         title(['L^p = ' num2str(p) ', L^q = ' num2str(q)]);
            %         drawnow;
        end
        fprintf('.');
        save([savefilename '_partial']);
    end
end
save(savefilename);
disp(['Saved as ' savefilename])
end


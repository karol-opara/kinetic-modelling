function runLowConcComparison(lambda,name)
%twoExperiments();
manyExperiments(lambda,name);
end

function manyExperiments(lambda,name)
load saveExperimentalData20130319_8of13experiments_NRTLvalidation

savefilename = ['Results/' 'save_LowConcComparison_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];

%lambda = [1, 0];

type = 'batch';

tic

p = {2,     'log',  1,  2}; % inner
q = {NaN,   NaN,    1,  'log'}; % outer loss (NaN means using only inner loss)

model = cell(8,6);
for i=1:length(p)
    parfor j = 1:8
        model{i,j} = EstimateKineticModel(data{j},p{i},q{i},type,lambda);
    end
end

save(savefilename);
disp(['Saved as ' savefilename]);
toc
end


function twoExperiments()
load saveExperimentalData20120807
data_OK = data{7};
data_err = data{5};

savefilename = ['Results/' 'save_LowConcComparison_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];

lambda = [1, 0];

type = 'batch';

tic

p = 2;
q = NaN;
control(1) = struct('p',p,'q',q);

p = 1;
q = NaN;
control(2) = struct('p',p,'q',q);

p = 0.5;
q = NaN;
control(3) = struct('p',p,'q',q);

p='log';
q = NaN;
control(4) = struct('p',p,'q',q);

p=2;
q = 'log';
control(5) = struct('p',p,'q',q);

p=1;
q = 'log';
control(6) = struct('p',p,'q',q);

m_ok=cell(6,1);
m_err=cell(6,1);

parfor i=1:6
    m_ok{i} = EstimateKineticModel(data_OK,control(i).p,control(i).q,type,lambda);
    m_err{i} = EstimateKineticModel(data_err,control(i).p,control(i).q,type,lambda);
end

m2_ok = m_ok{1};
m2_err = m_err{1};
m1_ok = m_ok{2};
m1_err = m_err{2};
m05_ok = m_ok{3};
m05_err = m_err{3};
mlog_ok = m_ok{4};
mlog_err = m_err{4};
mlog_2_ok = m_ok{5};
mlog_2_err = m_err{5};
mlog_1_ok = m_ok{6};
mlog_1_err = m_err{6};

save(savefilename);
disp(['Saved as ' savefilename]);
toc
end
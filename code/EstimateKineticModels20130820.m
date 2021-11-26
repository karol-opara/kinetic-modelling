function model = EstimateKineticModels20130820(name)
data=[];
load('saveExperimentalData20130820_batch');
savefilename = ['Results/' 'save_KineticModel_Batch_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];

weights = [0 1 1 1 1 1];

lambda=[0.01,0.03,0.05,0.08,0.1,0.3,0.7,1.5,2,0.005];
model=cell(length(data),length(lambda));
for id = 1:length(data)
    for i=1:length(lambda)
        Lambda=[1, lambda(i)];
        
        q = 'log'; % outer (inter-component) loss funciton
        p = 2; % inner loss function
        % for standard fit use
        %  p = 2;
        %  q = NaN;
        
        options.MaxFunEvals = 1e5;
        
        model{id,i} = EstimateKineticModel(data{id},p,q,'batch',Lambda,options,weights);
        
        fprintf(['Iteration ' num2str(i) ', lambda = ' num2str(Lambda(2))]);
        model{id,i}.k
        
        figure(i)
        plotKineticModelFit(model{id,i});
        axis([-1, max(data{id}.timeZ)+10, -0.2, Inf]);
        title(['Data ID = ' num2str(id) ', lambda = ' num2str(Lambda(2))]);
    end
end

save(savefilename);
end
function processEstimateKineticModels20130820
load('../data/save_KineticModel_Batch_2013-08-20_202623_allData_regularized_LLogL2')

for id = 1:length(data)
    for i=1:length(lambda)
        figure()
        plotKineticModelFit(model{id,i});
        axis([-1, max(data{id}.timeZ)+10, -0.2, Inf]);
        title(['Data ID = ' num2str(id) ', lambda = ' num2str(Lambda(2))]);
    end
end
end
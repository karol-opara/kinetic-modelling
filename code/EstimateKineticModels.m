function models = EstimateKineticModels(name)
%load saveExperimentalData20120726;
%load saveExperimentalData20120807;
load saveExperimentalData20130319_8of13experiments_NRTLvalidation

savefilename = ['Results/' 'save_KineticModels_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];
%p=[0.25 0.5 1 2];
p=0.5;
for i = 1:length(data)
    for j = 1:length(p)
        models{i,j}=EstimateKineticModel(data{i},p(j),0.5);
        save(savefilename);
    end
end

end
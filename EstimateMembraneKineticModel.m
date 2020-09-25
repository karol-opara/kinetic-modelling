function model = EstimateMembraneKineticModel(name,lambda,id)
dataRetentate=[];
load('saveExperimentalData20130819_membrane')

savefilename = ['Results/' 'save_KineticModel_Membrane_' ...
    datestr(now,'yyyy-mm-dd_HHMMSS') '_' name];
q = 'log'; % outer (inter-component) loss funciton
p = 2; % inner loss function
% for standard fit use
%  p = 2;
%  q = NaN;

dataRetentate{id}.z0ml(1) = 15;

options.MaxFunEvals = 1e3;
model=EstimateKineticModel(dataRetentate{id},p,q,'membrane',lambda,options);
model.k

plotKineticModelFit(model);
axis([-1, max(dataRetentate{id}.timeZ)+10, -0.2, Inf]);

save(savefilename);

end
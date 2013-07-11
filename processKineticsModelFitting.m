load Results\save_KineticModels_2012-08-16_172142_DataFrom20120807_L2_Lp
%load Results/save_KineticModels_2012-08-17_152529_DataFrom20120807_L2_Lp_studentized
%close all
[r c]=size(models);
for j=4%1:c
    %figure()
    for i =1:length(models)
        %kHandmade = getHandmadeKineticModelFit();
        %k=kHandmade{i};
        %k([1 3 5]) = k([1 3 5])/2;
        
        %subplot(2,2,i-3)
        figure()
        plotKineticModelFit(data{i}.timeZ,data{i}.zml, k, data{i}.z0ml,true,'-');
%         title(['L' num2str(p(j)) ' experiment ' num2str(i)...
%             ', NaOH = ' num2str(data{i}.NaOH)...
%             ', TG:MeOH = ' num2str(data{i}.TGMeOHRatio)]);
        
%axis([-Inf Inf 0 2.5]);
        models{i}.k
    end
end

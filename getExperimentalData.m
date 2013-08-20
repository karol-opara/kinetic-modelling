function [dataRetentate, dataPermeate]=getExperimentalData(type)
if (nargin==0)
    warning('getExperimentalData:NoType','No type provided, used batch as default');
    type = 'batch';
end
if(strcmp(type,'batch'))
    %dataRetentate = getBatchReactionData20130319();
    dataRetentate = getBatchReactionData20130820();
    dataPermeate=[];
elseif (strcmp(type,'membrane'))
    [dataPermeate, dataRetentate] = getMembraneReactionData();
else
    error('getExperimentalData:WrongType','Unsupported type')
end
end

function [dataPermeate, dataRetentate] = getMembraneReactionData()
expDataP = readMembraneDataPermeate();
expDataR = readMembraneDataRetentate();
%expData = expData([2,6,7,8,9,11,12,13]); % remove bad-looking data

wtpDataP = expData2wtp(expDataP);
wtpzDataP = calculateZ(wtpDataP);
mDataP = wtpData2mpml(wtpzDataP);
%dataP = removeWrongPoints(mDataP);
dataPermeate = mDataP;
close all
%plotWtpData(data);
%plotMpData(data);
plotMlData(dataPermeate, ' Permeate');




wtpDataR = expData2wtp(expDataR);
wtpzDataR = calculateZ(wtpDataR);
mDataR = wtpData2mpml(wtpzDataR);
%dataR = removeWrongPoints(mDataR);
dataRetentate = mDataP;
plotMlData(dataRetentate, ' retentate');


%save saveMembraneExperimentalData20130711


end

function expData = readMembraneDataPermeate()
table = xlsread('../Data calculation/2013-07-11_data_ppmwt_mr_permeate.xlsx');             %original data x,y in ppm,wt%
noOfExperiments = 4;
NaNLines = [1 9 15 25 35];
xOffset = 8;
for i = 1:noOfExperiments
    ex = table((1+NaNLines(i)):(NaNLines(i+1)-1),:);
    expData{i} = struct('timeX', ex(:,1),'timeY', ex(:,1),'timeZ', ex(:,1),...
        'temperature', ex(1,2), 'TGMeOHRatio', ex(1,3), 'NaOH', ex(1,4),...
        'y1wtp', ex(:,5), 'y26ppm', ex(:,6:10),...
        'yVolume', ex(:,11), 'yDensity', ex(:,12),...
        'x1wtp', ex(:,xOffset+5), 'x26ppm', ex(:,xOffset+(6:10)),...
        'xVolume', ex(:,xOffset+11), 'xDensity', ex(:,xOffset+12),...
        'z0ml', table(NaNLines(i),23:28),...
        'reactorVolume',ex(:,22)*1e3);
    expData{i}.zVolume = expData{i}.xVolume + expData{i}.yVolume;
end
end

function expData = readMembraneDataRetentate()
table = xlsread('../Data calculation/2013-07-11_data_ppmwt_mr_retentate.xlsx');             %original data x,y in ppm,wt%
noOfExperiments = 4;
NaNLines = [1 9 15 25 35];
xOffset = 8;
for i = 1:noOfExperiments
    ex = table((1+NaNLines(i)):(NaNLines(i+1)-1),:);
    expData{i} = struct('timeX', ex(:,1),'timeY', ex(:,1),'timeZ', ex(:,1),...
        'temperature', ex(1,2), 'TGMeOHRatio', ex(1,3), 'NaOH', ex(1,4),...
        'y1wtp', ex(:,5), 'y26ppm', ex(:,6:10),...
        'yVolume', ex(:,11), 'yDensity', ex(:,12),...
        'x1wtp', ex(:,xOffset+5), 'x26ppm', ex(:,xOffset+(6:10)),...
        'xVolume', ex(:,xOffset+11), 'xDensity', ex(:,xOffset+12),...
        'z0ml', table(NaNLines(i),23:28),...
        'reactorVolume',ex(:,22)*1e3);
    expData{i}.zVolume = expData{i}.xVolume + expData{i}.yVolume;
end
end


function data = getBatchReactionData20130319()
expData = readData20130319();

expData = expData([2,6,7,8,9,11,12,13]); % remove bad-looking data

wtpData = expData2wtp(expData);
wtpzData = calculateZ(wtpData);
mData = wtpData2mpml(wtpzData);

data = mData;
data = removeWrongPoints(mData);

%save saveExperimentalData20130319_8of13experiments_NRTLvalidation

%close all
%plotWtpData(data);
%plotMpData(data);
plotMlData(data);
end


function data = getBatchReactionData20130820()
expData = readData20130820();
wtpData = expData2wtp(expData);
wtpzData = calculateZ(wtpData);
mData = wtpData2mpml(wtpzData);
mData{1} = removeSinglePoint(mData{1}, 6, 'x');
mData{1} = removeSinglePoint(mData{1}, 6, 'y');
mData{1} = removeSinglePoint(mData{1}, 6, 'z');
data = mData;
plotMlData(data);

%save('saveExperimentalData20130820_batch','data');
end

function expData = readData20130820()
table = xlsread('../Data calculation/2013-08-06_data_ppmwt.xlsx');             %original data x,y in ppm,wt%
noOfExperiments = 7;
NaNLines = [1 11 23 30 39 48 59 70];
xOffset = 8;
for i = 1:noOfExperiments
    ex = table((1+NaNLines(i)):(NaNLines(i+1)-1),:);
    expData{i} = struct('timeX', ex(:,1),'timeY', ex(:,1),'timeZ', ex(:,1),...
        'temperature', ex(1,2), 'TGMeOHRatio', ex(1,3), 'NaOH', ex(1,4),...
        'y1wtp', ex(:,5), 'y26ppm', ex(:,6:10),...
        'yVolume', ex(:,11), 'yDensity', ex(:,12),...
        'x1wtp', ex(:,xOffset+5), 'x26ppm', ex(:,xOffset+(6:10)),...
        'xVolume', ex(:,xOffset+11), 'xDensity', ex(:,xOffset+12),...
        'z0ml', table(NaNLines(i),23:28),...
        'reactorVolume',ex(:,22)*1e3);
    expData{i}.zVolume = expData{i}.xVolume + expData{i}.yVolume;
end
end

function expData = readData20130319()
table = xlsread('../Data calculation/2013-01-24_data_ppmwt.xlsx');             %original data x,y in ppm,wt%
noOfExperiments = 13;
NaNLines = [1 14 27 37 48 60 72 84 96 103 112 121 132 143];
xOffset = 8;
for i = 1:noOfExperiments
    ex = table((1+NaNLines(i)):(NaNLines(i+1)-1),:);
    expData{i} = struct('timeX', ex(:,1),'timeY', ex(:,1),'timeZ', ex(:,1),...
        'temperature', ex(1,2), 'TGMeOHRatio', ex(1,3), 'NaOH', ex(1,4),...
        'y1wtp', ex(:,5), 'y26ppm', ex(:,6:10),...
        'yVolume', ex(:,11), 'yDensity', ex(:,12),...
        'x1wtp', ex(:,xOffset+5), 'x26ppm', ex(:,xOffset+(6:10)),...
        'xVolume', ex(:,xOffset+11), 'xDensity', ex(:,xOffset+12),...
        'z0ml', table(NaNLines(i),23:28),...
        'reactorVolume',ex(:,22)*1e3);
    expData{i}.zVolume = expData{i}.xVolume + expData{i}.yVolume;
end
end

function wtpData = expData2wtp(expData)
for i = 1:length(expData)
    wtpData{i} = rmfield(expData{i}, {'y1wtp', 'x1wtp', 'y26ppm', 'x26ppm'});
    sumxppm = repmat(sum(expData{i}.x26ppm.').',1,5);
    xratio = expData{i}.x26ppm./sumxppm;
    wtpData{i}.xwtp = [expData{i}.x1wtp repmat((100-expData{i}.x1wtp),1,5).*xratio];
    
    sumyppm = repmat(sum(expData{i}.y26ppm.').',1,5);
    yratio = expData{i}.y26ppm./sumyppm;
    wtpData{i}.ywtp = [expData{i}.y1wtp repmat((100-expData{i}.y1wtp),1,5).*yratio];
end
end

function mData = wtpData2mpml(wtpData)
for i = 1:length(wtpData)
    data = wtpData{i};
    wtData{i} = data;
    
    [~, c] = size(data.xwtp);
    
    xmass = repmat(data.xDensity.*data.xVolume,1,c);
    ymass = repmat(data.yDensity.*data.yVolume,1,c);
    
    xwt = 1e-2*data.xwtp.*xmass;
    ywt = 1e-2*data.ywtp.*ymass; % weights of each component in the experimetal sample
    xmole = wt2m(xwt);
    ymole = wt2m(ywt);
    sumxmole = repmat(sum(xmole.').',1,c);
    sumymole = repmat(sum(ymole.').',1,c);
    
    data.xmp = 100*xmole./sumxmole;
    data.ymp = 100*ymole./sumymole;
    
    data.xml = xmole./(1e-3*repmat(data.xVolume,1,c)); % mole per litre
    data.yml = ymole./(1e-3*repmat(data.yVolume,1,c));
    
    data.zwtBulk = (xwt+ywt).*repmat(data.reactorVolume./data.zVolume,1,6);
    
    mData{i} = data;
end
end



function wtpzData = calculateZ(wtpData)
for i = 1:length(wtpData)
    data = wtpData{i};
    wtpzData{i} = data;
    
    [~, c] = size(data.xwtp);
    
    xmass = repmat(data.xDensity.*data.xVolume,1,c);
    ymass = repmat(data.yDensity.*data.yVolume,1,c);
    
    zwt = 1e-2*(data.xwtp.*xmass+data.ywtp.*ymass); % weights of each component in the experimetal sample
    zmole = wt2m(zwt); % convert weights to moles
    zVolume = 1e-3*(data.xVolume+data.yVolume); % in litres
    zmoleL = zmole./repmat(zVolume,1,c);
    wtpzData{i}.zwtp = 1e2*zwt./(xmass+ymass);
    wtpzData{i}.zmp = wtp2mp(wtpzData{i}.zwtp);
    wtpzData{i}.zml = zmoleL;
end
end

% function mpData = wtpData2mp(wtpzData)
% for i = 1:length(wtpzData)
%     mpData{i} = wtpzData{i};
%     
%     mpData{i}.xmp = wtp2mp(wtpzData{i}.xwtp);
%     mpData{i}.ymp = wtp2mp(wtpzData{i}.ywtp);
%     mpData{i}.zmp = wtp2mp(wtpzData{i}.zwtp);
% end
% end

function [m] = wt2m(wt)
% convert weights to moles
% Molecural weights (sometimes averaged)
MW=[32      % MeOH
    92.1    % GLY
    849.5   % TG
    597     % DG
    344.5   % MG
    284.5   % FAME
    ];
[r, ~] = size(wt);
m = wt./repmat(MW.',r,1); % number of moles
end

function [mp ml] = wtp2mp(wtp)
% Convert weight percent to mole percent
    [~, c] = size(wtp);
    unitMass = 1;
    mole = wt2m(wtp*unitMass);
    summole = repmat(sum(mole.').',1,c);
    mp = mole./summole*100;
end

function plotWtpData(wtpData)
for i = 1:length(wtpData)
figure()
[hx, hy, hz] = plotConcentrations(wtpData{i}.timeX, wtpData{i}.xwtp,...
    wtpData{i}.timeY, wtpData{i}.ywtp, wtpData{i}.timeZ, wtpData{i}.zwtp);
xlabel(hx,'time [min]');
xlabel(hy,'time [min]');
xlabel(hz,'time [min]');
ylabel(hx,'x [wt%]');
ylabel(hy,'y [wt%]');
ylabel(hz,'z [wt%]');
title(hy,['experiment ' num2str(i)...
    ', NaOH = ' num2str(wtpData{i}.NaOH)...
    ', TG:MeOH = ' num2str(wtpData{i}.TGMeOHRatio)]);
end
end

function plotMpData(mpData)
for i = 1:length(mpData)
figure()
[hx, hy, hz] = plotConcentrations(mpData{i}.timeX, mpData{i}.xmp,...
    mpData{i}.timeY, mpData{i}.ymp, ...
    mpData{i}.timeZ, mpData{i}.zmp);
xlabel(hx,'time [min]');
xlabel(hy,'time [min]');
xlabel(hz,'time [min]');
ylabel(hx,'x [mole%]');
ylabel(hy,'y [mole%]');
ylabel(hz,'z [mole%]');
title(hy,['experiment ' num2str(i)...
    ', NaOH = ' num2str(mpData{i}.NaOH)...
    ', TG:MeOH = ' num2str(mpData{i}.TGMeOHRatio)]);
end
end

function plotMlData(mlData, name)
if (nargin == 1)
    name = '';
end
for i = 1:length(mlData)
figure()
[hx, hy, hz] = plotConcentrations(mlData{i}.timeX, mlData{i}.xml,...
    mlData{i}.timeY, mlData{i}.yml, ...
    [0; mlData{i}.timeZ], [mlData{i}.z0ml; mlData{i}.zml]);
mlData{i}.z0ml
xlabel(hx,'time [min]');
xlabel(hy,'time [min]');
xlabel(hz,'time [min]');
ylabel(hx,'x [mole/l]');
ylabel(hy,'y [mole/l]');
ylabel(hz,'z [mole/l]');
title(hy,['experiment ' num2str(i)...
    ', NaOH = ' num2str(mlData{i}.NaOH)...
    ', TG:MeOH = ' num2str(mlData{i}.TGMeOHRatio)...
    name]);
end
end

function data = removeWrongPoints(data)
data{4} = removeSinglePoint(data{4}, 6, 'x');
data{4} = removeSinglePoint(data{4}, 6, 'y');
data{4} = removeSinglePoint(data{4}, 6, 'z');

data{4} = removeSinglePoint(data{4}, 7, 'x');
data{4} = removeSinglePoint(data{4}, 7, 'y');
data{4} = removeSinglePoint(data{4}, 7, 'z');

data{3} = removeSinglePoint(data{3}, 6, 'x');
data{3} = removeSinglePoint(data{3}, 6, 'z');
data{3} = removeSinglePoint(data{3}, 7, 'y');
data{3} = removeSinglePoint(data{3}, 7, 'z');
end

function data = removeSinglePoint(data, id, type)
if (type == 'x')
    data.timeX(id)=[];
    data.xVolume(id) = [];
    data.xDensity(id) = [];
    data.xwtp(id,:) = [];
    data.xmp(id,:) = [];
    data.xml(id,:) = [];
elseif (type == 'y')
    data.timeY(id)=[];
    data.yVolume(id) = [];
    data.yDensity(id) = [];
    data.ywtp(id,:) = [];
    data.ymp(id,:) = [];
    data.yml(id,:) = [];
elseif (type == 'z')
    data.timeZ(id)=[];
    data.zwtp(id,:) = [];
    data.zmp(id,:) = [];
    data.zml(id,:) = [];
else
    warning('getExperimentalData:removeSinglePoint','Wrong point type specified');
end
end
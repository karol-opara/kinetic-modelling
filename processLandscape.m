function processLandscape(path,plotType)
load(path);

[r, c] = size(K);
for i = 1:r
    for j = 1:c
        k = KProjected{i,j};
        kF(i,j) = sum(k([1 3 5]));
        kB(i,j) = sum(k([2 4 6]));
    end
end


% get coordinates of "the best point"
[errMinCols, ri] = min(ErrLp);
[~, idxBckBest] = min(errMinCols);
idxFwdBest = ri(idxBckBest);

if (strcmp(plotType,'contourf'))
    contourf(kF,kB,log(1+ErrLp),18);
    mean(mean(ErrLp))
    hl = line(kF,kB);
    set(hl,'LineStyle','none','Marker','.','Color','k');
    xlabel('forward: k_1 + k_3 + k_5')
    ylabel('backward: k_2 + k_4 + k_6')
    set(gca,'XScale','log','YScale','log');
    % title(['TG:MeOH = 1:' num2str(mData{condition}.TGMeOHRatio)...
    %     ', NaOH = ' num2str(mData{condition}.NaOH) ', Norm L' num2str(p)]);
    title(['Norm L' num2str(p)]);
    
    K = getLitertureK();
    h = plotK(K,[1 1 1]);
    
    hl = line(kF(idxFwdBest,idxBckBest), kB(idxFwdBest,idxBckBest));
    set(hl,'LineStyle','none','Marker','d','Color',[1 1 1]);
elseif(strcmp(plotType,'plot'))
    zTime = data.timeZ;
    zml = data.zml;
    
    kij = KProjected{idxFwdBest,idxBckBest}
    plotKineticModelFit(zTime,zml,KProjected{idxFwdBest,idxBckBest},data.z0ml,'batch');
%     title(['k_f = ' num2str(sum(kij([1 3 5]))) ', k_b = ' num2str(sum(kij([2 4 6])))...
%         ' Norm L' num2str(p)]);
    title([' Norm L' num2str(p)]);
    %ViewPlot.Save(no,['KineticModelFit' num2str(no)]);
    if (zml(1,1)==13.058476552736959) % just a hack for preparing manuscript
        axis([-1 70 -0.03 1.5]);
    else
        axis([-1 70 -0.05 2.5]);
    end
else
    warning('processLandscape:main','Unsupported plot type');
end
end

function hl = plotK(K,color)
kFLit = zeros(1,length(K));
kBLit = zeros(1,length(K));
labels = {};
for i=1:length(K)
    kFLit(i) = K{i}.kf;
    kBLit(i) = K{i}.kb;
    labels{i}=K{i}.id;
end

hl = line(kFLit,kBLit);
set(hl,'Marker','o','LineStyle','none','Color',color);
text(kFLit, kBLit, labels, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right','Color',color);
end

function plotSingleResult(zTime, zml, K, Z0, kF, kB, ErrLp, i,j)
plotKineticModelFit(zTime, zml, K{i,j}, Z0{i,j});
title(['Kinetic model fit; k_f = ' num2str(kF(i,j),'%1.2f')...
    ', k_b = ' num2str(kB(i,j),'%1.2f') ', err = ' num2str(ErrLp(i,j),'%1.2f')]);
end

function K = getLitertureK()
i=1;
k = [0.0500
    0.1100
    0.2150
    1.2280
    0.2420
    0.0070];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','1');
i=i+1;

k=[0.1030
    0.0310
    0.0630
    0.0100
    0.0160
    0.1750];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','2');
i=i+1;

k=[0.0286
    0.0144
    0.0058
    0.0213
    0.0111
    0.0005
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','3');
i=i+1;

k=[0.8000
    5.9500
    10.5000
    15.9000
    0.3400
    0.0035];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4a');
i=i+1;

k=[1.5500
    8.5000
    20.5000
    22.5000
    0.6100
    0.0012];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4b');
i=i+1;

k=[2.0500
    10.9000
    30.1000
    29.5000
    0.8300
    0.0001];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4c');
i=i+1;

k=[1.5000
    13.7000
    23.0000
    41.4000
    0.4000
    0.0026
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4d');
i=i+1;

k=[3.0600
    23.9000
    32.5000
    57.5000
    0.5400
    0.0009
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4e');
i=i+1;

k=[4.0000
    27.0000
    55.0000
    65.5000
    0.9100
    0.0001
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','4f');
i=i+1;

k=[2.5790
    0.0200
    0.6000
    0.1010
    0.9000
    0.0210];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','5a');
i=i+1;

k=[2.6000
    0.2480
    1.1860
    0.2270
    2.3030
    0.0220
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','5b');
i=i+1;

k=[2.6200
    0.7000
    1.2100
    0.4000
    2.3600
    0.0280
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','5c');
i=i+1;

k=[0.0247
    0.0558
    0.0703
    0.0015
    0.0394
    0.0042];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','6a');
i=i+1;
k=[0.0772
    0.1680
    0.0972
    0.0265
    0.0670
    0.0088];

K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','6b');
i=i+1;
k=[0.0443
    0.2334
    0.0645
    0.0699
    0.2681
    0.0047];

K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','6c');
i=i+1;

k=[0.0879
    0.4777
    0.1555
    0.1396
    0.7478
    0.0061];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','6d');
i = i+1;


k=[0.01057
    0.00000
    0.11840
    0.08187
    0.13100
    0.00201
    ];
K{i} = struct('k',k,'kf',sum(k([1 3 5])), 'kb',sum(k([2 4 6])),'id','7');

end
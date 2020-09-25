err = [1.8, 0.22, 0.64, 0.19, 0.03, 0.01];
x = 0:0.01:2;
L2 = @(x) x.^2;
L1 = @(x) abs(x);
L05 = @(x) sqrt(abs(x));
LLog = @(x) log1p(abs(x));

subplot(2,2,[1 3])
h=plot(x,L2(x),'-',x,L1(x),'--',x,L05(x),'-.',x,LLog(x),':');
set(h,'LineWidth',1.2);
set(gca,'XGrid','on','YGrid','on');
l = legend('Square','Absolute','Root','Log','Location','NW');
%l = legend('Square $e^2$','Absolute $e$','Root $\sqrt{e}$','Log $\log(1+e)$','Location','NW');
%set(l,'Interpreter','Latex');
axis([0 2 0 4]);
title('a) loss functions')
xlabel('error');
ylabel('loss');

subplot(2,2,2)
bar(diag(err),'stacked');
set(gca,'YGrid','on','XLim',[0 7]);
title('b) example errors')
ylabel('error')
xlabel('observation number')

subplot(2,2,4)
barh(100*[ LLog(err)/sum(LLog(err)); L05(err)/sum(L05(err)); L1(err)/sum(L1(err)); L2(err)/sum(L2(err))],'stacked')
axis([0 100 -Inf Inf]);
title('c) relative influence')
xlabel('influence [%]')
set(gca,'yticklabel',{'Log','Root','Absolute','Square'});
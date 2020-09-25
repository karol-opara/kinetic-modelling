x=-20:0.5:20;
[X,Y]=meshgrid(x,x);
Z = log1p(abs(X))+0.3*log1p(abs(Y));
subplot(2,1,1)
surfc(X,Y,Z)
ylabel('Regularization term')
xlabel('Fit term')
zlabel('log(1+fit)+0.3*log(1+reg)')
subplot(2,1,2)
surfc(X,Y,Z)
ylabel('Regularization term')
xlabel('Fit term')
zlabel('log(1+fit)+0.3*log(1+reg)')
view([-30 80]);
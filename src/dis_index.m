clear
clc
close all
global lg
dnorm = load('D:\Paper 2\te_process\d00.dat');   dnorm = dnorm(1:480,lg); 
d01 = load('D:\Paper 2\te_process\d15.dat');     d01 = d01(:,lg);
[dnorm,dnormx] = mncn(dnorm);
[d01,d01x] = mncn(d01);
n1 = size(dnorm,1); n2 = size(d01,1);
R = n1/(n1+n2)*cov(dnorm)+n2/(n1+n2)*cov(d01);
[v,s] = eig(R);
y1 = sqrt(n1/(n1+n2))*dnorm*v*s^(-0.5);
[yv1,ys1] = eig(cov(y1));
D1 = 4/(size(dnorm,2))*sum((diag(ys1)-0.5).^2);

d01 = load('D:\Paper 2\te_process\d15_te.dat');   d01 = d01(:,lg);
Dte = scale(d01,dnormx(1,:),dnormx(2,:));
for i = 1:size(Dte,1)
   Samp = Dte(i,:);
   Dte_y = Samp*v*s^(-0.5);
   [~,ys1] = eig(cov(Dte_y));
   D2(i) = 4/(size(Dte,2))*sum((diag(ys1)-0.5).^2);
end
plot(D2)
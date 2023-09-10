close all
clear
clc
alpha = 2.43;
beta = 1;
pdf = @(x)gampdf(x,alpha,beta);
proppdf = @(x,y)gampdf(x,floor(alpha),floor(alpha)/alpha);
proprnd = @(x)sum(...
              exprnd(floor(alpha)/alpha,floor(alpha),1));
nsamples = 5000;
smpl = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,...
                'proppdf',proppdf);
[f,xi,w] = ksdensity(smpl,'support','positive','npoints',1000);
figure()
f2 = pdf(xi);
hold on
plot(xi,f2,'r')
plot(xi,f,'b')
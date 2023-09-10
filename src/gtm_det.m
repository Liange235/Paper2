clear all
close all
clc
d0_norm = load('D:\Paper 2\te_process\d00.dat'); 
data = d0_norm';
[D,Dx] = mncn(data);
% n = 20;
% cv = cvpartition(size(D,1), 'kfold', n);
% alphas = logspace(1,2,10);
% logprobs = zeros(10,numel(alphas));
% for i = 1:numel(alphas)
%   for j = 1:n
%     net = gtm_make(D(cv.training(j),:),'small', 'rbfgrid', [6 7], 'regul', alphas);
%     logprobs(j,i) = -sum(log(gtmprob(net, D(cv.test(j),:))));
%   end
% end
%% plot cross-validation results
% figure(1)
% semilogx(alphas,mean(logprobs),'bo-');
% the plot shows that alpha = XXX is en elbow point, thus suitable as
% regularization parameter

% [~,ind] = min(mean(logprobs));
% alpha = alphas(ind);
%% Use wrapper function gtm_make to train a GTM
net = gtm_make(D,'normal', 'rbfgrid', [6 7],'regul', 0);
d01 = load('D:\Paper 2\te_process\d05_te.dat'); 
Dte = scale(d01,Dx(1,:),Dx(2,:));
[post,a,b] = gtmpost(net, D);  [post_te,a_te,b_te] = gtmpost(net,Dte); a = log(a); a_te = log(a_te);
post_norm = mean(a)
% [post,a,b] = gtmpost(net, Dte);
sub = a_te-ones(size(post_te,1),1)*post_norm;
Dism = diag(sqrt(sub*sub'));
figure()
plot(Dism,'r');


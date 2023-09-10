clear all
close all
clc
norm = load('D:\yanyixia\TE_process\d00.dat'); norm = norm(floor(linspace(size(norm,1),1,50)),:); label_0 = 1*ones(size(norm,1),1);
fault1 = load('D:\yanyixia\TE_process\d04.dat'); fault1 = fault1(floor(linspace(size(fault1,1),1,50)),:); label_1 = 2*ones(size(fault1,1),1);
fault2 = load('D:\yanyixia\TE_process\d05.dat'); fault2 = fault2(floor(linspace(size(fault2,1),1,50)),:); label_2 = 3*ones(size(fault2,1),1);
data = [norm;fault1;fault2]; labels = [label_0;label_1;label_2];
sD = som_normalize(data, 'var');
cv = cvpartition(size(sD,1), 'kfold', 10);
alphas = logspace(-1,-4,10);
logprobs = zeros(10,numel(alphas));
for i = 1:numel(alphas)
  for j = 1:10
    net = gtm_make(sD, 'rbfgrid', [3 3], 'regul', alphas(i));
    logprobs(j,i) = -sum(log(gtmprob(net, sD)));
  end
end
%% plot cross-validation results
figure(1)
semilogx(alphas,mean(logprobs),'bo-');
% the plot shows that alpha = XXX is en elbow point, thus suitable as
% regularization parameter
[~,ind] = min(mean(logprobs));
alpha = alphas(ind);
%% Use wrapper function gtm_make to train a GTM
net = gtm_make(sD, 'rbfgrid', [3 3],'regul', alpha);
gtm_show(net, 'mags', 'data', sD, 'groups', labels);
gtm_show(net, 'mags', 'data', sD, 'groups', labels,'mapping','mode');

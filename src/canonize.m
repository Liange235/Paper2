function [coef,S_cat] = canonize(X,grp) 
%********************CVA 规范变量分析方法***********************************
% Performs canonical variate matrix of data in X. 
% gcb, 04 March 2006 
% X: N x D data matrix, data points in rows, variables in columns 
% grp: grouping vector, containing integers between 1 and K 
% where K is the number of groups 
% coef: coefficients for forming canonical variates 
% (eigenvectors of inv(S)*B) 
% score: canonical variate scores 
% S: within-groups covariance matrix 
% B: between-groups covariance matrix 
%**************************************************************************
N = size(X,1); 
% number of data (rows) 
D = size(X,2); 
% number of variables (columns) 
K = max(grp); % number of groups 


xmg = mean(X); % global mean vector (1 x D) 
xmk = zeros(K,D); % will hold group means 
nk = zeros(K,1); % number of data per group 
for k = 1:K 
xmk(k,:) = mean(X(grp==k,:)); 
nk(k) = size(X(grp==k,:),1); 
end % calc within-groups cov. mat. 
S = zeros(D,D);
for i = 1:K
    xdiff = X(grp==i,:)-ones(size(X(grp==i,:),1),1)*xmk(i,:);
    S_cat{i} = xdiff'*xdiff;
    S = S_cat{i}+S;
end
S = S/(N-K); 
% calc between-groups cov. mat.
B = zeros(D,D); 
for k = 1:K 
B = B + nk(k)*(xmk(k,:)-xmg)'*(xmk(k,:)-xmg); 
end 
B = (K/(K-1))*B/N; 
% do eigenanalyis and put in descending order 
[coef,ev] = eig(B,S); 
[ev,iev] = sort(diag(ev),1,'descend'); 
coef = coef(:,iev); 





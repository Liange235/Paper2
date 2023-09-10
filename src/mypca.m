%% PCA 去相关性&降维
% Calculate the covariance matrix
function Wxpca = mypca(D)
if size(D,1) >= size(D,2)
    C = D' * D;	 
else
    C = D  * D';	 
end
[eigVecs, eigVals] = eig(C);
eigVals = abs(diag(eigVals));
[~,index] = sort(eigVals,'descend');
eigVals = eigVals(index);
eigVecs = eigVecs(:,index);
sumEigVal = sum(eigVals);
for i = 1:size(eigVals,1)
   if (sum(eigVals(1:i)/sumEigVal))>= 0.95;
       ind = i;
       break
   end
end
eigVals(ind+1:end) = [];
eigVecs(:,ind+1:end) = [];
if size(D,1) >= size(D,2)
    Wxpca = eigVecs;
else
    Wxpca = D' * eigVecs * diag(1 ./ sqrt(eigVals));
end

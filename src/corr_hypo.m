%% 采用传统线性相关系数以及p值检验相关显著性程度
function [y,y_alt] = corr_hypo(x,shift)
global mar
[~,pval] = corr(x);
pval = circshift(pval,-shift); % 向下循环移动1位
n = size(x,2);
y_all = [1:size(x,2)];
ind = 1;
y_alt = 1;
a = cell(1,n);
for i = 1:n
    a{i} = find(pval(i,i+1:end)<=0.05)+i;
end
while(1)  
   y_tem = [y_alt,a{ind}];
   y_all = y_all(~ismember(y_all ,y_tem));
    if (isempty(y_all))
       break;
   end
   ind = y_all(1);
   y_alt = [y_alt,ind];
end
%%%%%%%%%%%%********变量数化成统一mar**********%%%%%%%%%%%%%%
i = 1;
while (numel(y_alt)<mar)
   y_tem = a{i};
   y_tem = y_tem(~ismember(a{i},y_alt));
   if isempty(y_tem)
       i = i+1;
       continue
   end
   y_alt = [y_alt,y_tem(1)];
   i = i+1;
end
y_alt = sort(y_alt);
y = x(:,y_alt);
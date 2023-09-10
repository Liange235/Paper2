function y = dynamic(x,t)
% 输入x和时延窗口t，对矩阵x，进行时延动态化
% 输出y
% 假设x belongs to mxn
% 1.首先设定统一的时延区间t对x进行动态化生成y，y belongs to (m-t)x[n(t+1)]
%  
[m,n] = size(x);
y  = [];
temp = 1:m-t;
lab = 1:m-t;
for i = 1:t
     temp = temp+ones(1,m-t);
     lab = [lab;temp];
end
lab  = lab';
for i = 1:n
    temp = x(:,i);
    y = [y temp(lab)];
end
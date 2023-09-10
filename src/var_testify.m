clear
close all
clc
global mar 
%********training*****************%
D = []; tao = 0;dlim = 480;  
lib = [0:21]'; numl = length(lib);
for i = 1:numl
str = num2str(lib(i));
if (lib(i)<10)
    str = strcat('0',num2str(lib(i)));
end
d_com = load(strcat('D:\Paper 2\te_process\d',str,'.dat'));
d_com = d_com(1:dlim,:);

   for j = 1:size(d_com,2)
      index = corr_hypo(d_com,j-1);
      var_n = size(index,2);
      Map(j,i) = var_n;
   end
end

%% %%*********for TE*******%%%%%%
% lib = sort(unique(Map(:,1)));lib_lg = zeros(1,numel(lib));
% for j = 1:numel(lib)
%   lg = 0;
%   for i = 1:22
%     temp = Map(:,i);
%       if (sum(temp==lib(j)))
%          lg = lg+1;
%       end
%   end
%   lib_lg(j) = lg;
% end
% plot(lib_lg,'r*-')
% set(gca,'XTickLabel',num2str(lib));

%% %%*********for other data*******%%%%%%
mar = max(min(Map));
global shift
for i = 1:size(Map,2)
    shift(i) = min(find(Map(:,i)==mar))-1;
end
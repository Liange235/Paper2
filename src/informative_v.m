clear
clc
close all
global lg
dnorm = load('D:\Paper 2\te_process\d00.dat'); D = dnorm;
[D,Dx] = mncn(D);
d01 = load('D:\Paper 2\te_process\d15_te.dat'); 
d01 = scale(d01,Dx(1,:),Dx(2,:));
Dte1 = d01(1:160,:);  Dte2 = d01(161:end,:);
a = 2;
for i = 1:52
     J(i) = mean(log(cosh(a*Dte1(:,i)))/a)-mean(log(cosh(a*Dte2(:,i)))/a);
%    J(i) = mean(exp(-0.99*Dte1(:,i).^2/2))-mean(exp(-0.99*Dte2(:,i).^2/2));
%    J(i) = mean(Dte1(:,i).^4)-mean(Dte2(:,i).^4);
   J(i) = J(i)^2;
end
% J = (J-mean(J))/std(J);
figure(); plot(J)
[~,ind] = sort(J,'descend');
lg = ind(1:3);
figure(); hold on
plot(d01(:,lg)); legend(strcat('variable',num2str(lg(1))),strcat('variable',num2str(lg(2))),strcat('variable',num2str(lg(3))),strcat('variable',num2str(lg(4))));
up_lim = max(max(d01(:,lg))); floor_lim = min(min(d01(:,lg)));
plot([160,160],[floor_lim,up_lim],'--k','LineWidth',2);


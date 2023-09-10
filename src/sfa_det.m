clear
clc
close all
dnorm = load('E:\My_Education\Graduate\Paper 2\te_process\d00.dat');    
[dnorm,Dx] = mncn(dnorm);  
d01 = load('E:\My_Education\Graduate\Paper 2\te_process\d05.dat');  
D = scale(d01,Dx(1,:),Dx(2,:));
D = [dnorm;D];
% create a SFA object
n = 52;
hdl = sfa1_create(n);
% perform the preprocessing step
sfa_step(hdl, dnorm, 'preprocessing');
% sfa_step(hdl, dnorm, 'expansion');
% close the algorithm
sfa_step(hdl, [], 'sfa');

% compute the output signal
Dte_sfa = sfa_execute(hdl, D);



Dte1 = Dte_sfa(1:500,:);  Dte2 = Dte_sfa(501:end,:);
a = 1.5;
for i = 1:n
%      J(i) = mean(log(cosh(a*Dte1(:,i)))/a)-mean(log(cosh(a*Dte2(:,i)))/a);
   J(i) = mean(exp(-0.99*Dte1(:,i).^2/2))-mean(exp(-0.99*Dte2(:,i).^2/2));
%    J(i) = mean(Dte1(:,i).^4)-mean(Dte2(:,i).^4);
   J(i) = J(i)^2;
end
% J = (J-mean(J))/std(J);
figure(); plot(J)
[~,ind] = sort(J,'descend');
lg = ind(1:3);
figure(); hold on
plot(Dte_sfa(:,lg(1:3))); legend(strcat('variable',num2str(lg(1))),strcat('variable',num2str(lg(2))),strcat('variable',num2str(lg(3))));
up_lim = max(max(Dte_sfa(:,lg))); floor_lim = min(min(Dte_sfa(:,lg)));
plot([500,500],[floor_lim,up_lim],'--k','LineWidth',2);
title(num2str(J(ind(1))));
norm_sfa = sfa_execute(hdl, dnorm); 
% plot(norm_sfa(:,lg),'r'); legend('f','p','n')
%%%%*******使用SOM对故障3进行监控*******%%%%%
%% 对sfa特征提取之后的正常样本进行SOM建模
norm_sfa = norm_sfa(:,lg);
sDiris = som_data_struct(norm_sfa,'name','Te (train)');
sMap = som_make(sDiris,'small','tracking',1); 
h0 = som_hits(sMap,sDiris.data);
[~,ind] = sort(h0,'descend');
N_Center = sMap.codebook(ind(1),:);
bmus = som_bmus(sMap,sDiris,1);
sub = sDiris.data-ones(size(bmus,1),1)*N_Center;
Disim = diag(sqrt(sub*sub'));
for i = 2:28
    N_Center(i,:) = sMap.codebook(ind(i),:);
    sub = sDiris.data-ones(size(bmus,1),1)*N_Center(i,:);
    Disim = Disim+diag(sqrt(sub*sub'));
end
[fun,varx] = ksdensity(Disim,'bandwidth',0.2,'function','cdf');
control_lim = varx(min(find(fun>0.95)));
clear sDiris
%% 利用在正常样本下所建的SOM模型对故障5进行监控
d01 = load('E:\My_Education\Graduate\Paper 2\te_process\d05_te.dat'); 
d01 = scale(d01,Dx(1,:),Dx(2,:));
Dte_sfa = sfa_execute(hdl, d01);
Dte_sfa = Dte_sfa(:,lg);
sDiris = som_data_struct(Dte_sfa,'name','Te (test)');
for i = 1:size(Dte_sfa,1)
%       sub = sDiris.data(i,:)*Vx{bmus(i)}-norm_center;
%       Disim_te(i) = norm(sub);
%       Disim_te(i) = norm(sDiris.data(i,:)-norm_center);
        sub = ones(size(N_Center,1),1)*sDiris.data(i,:)-N_Center;
        Disim_te(i) = sum(diag(sqrt(sub*sub')));
end
figure(); hold on
plot([1:length(Disim_te)],Disim_te,'b');
plot([1:length(Disim_te)],control_lim,'r.-');
fdr = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; fdr = roundn(fdr,-2);
far = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  far = roundn(far,-2);
title(strcat('FAR:',num2str(far),'%','  FDR:',num2str(fdr),'%'),'FontName','Times','FontSize',13.5);
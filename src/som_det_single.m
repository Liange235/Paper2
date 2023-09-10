clear
clc
close all
dnorm = load('D:\Paper 2\te_process\d00.dat'); 
wid = 10;%*% window width
for i = 1:size(dnorm,1)-wid+1  D(i,:) = moving_window(dnorm,i+wid-1,wid); end
[D,Dx] = mncn(D);
sDiris = som_data_struct(D,'name','Te (train)');
sMap = som_make(sDiris,'small','tracking',1); 
h0 = som_hits(sMap,sDiris.data);
[~,ind] = sort(h0,'descend');
N_Center = sMap.codebook(ind(1),:);
% N_Center = mean(sMap.codebook(ind(1:2),:));
bmus = som_bmus(sMap,sDiris.data,1);
sub = sDiris.data-ones(size(bmus,1),1)*N_Center;
Disim = diag(sqrt(sub*sub'));
% for i = 2:10
%     N_Center(i,:) = sMap.codebook(ind(i),:);
%     sub = sDiris.data-ones(size(bmus,1),1)*N_Center(i,:);
%     Disim = Disim+diag(sqrt(sub*sub'));
% end
[fun,varx] = ksdensity(Disim,'bandwidth',0.1,'function','cdf');
control_lim = varx(min(find(fun>0.999)));
% control_lim = max(varx)+7;
clear sDiris

%%%***************************%%%%%
d01 = load('D:\Paper 2\te_process\d21_te.dat');
for i = 1:size(d01,1)-wid+1  Dte(i,:) = moving_window(d01,i+wid-1,wid); end
Dte = scale(Dte,Dx(1,:),Dx(2,:));
sDiris = som_data_struct(Dte,'name','Te (test)');
for j = 1:size(Dte,1)
    sub = ones(size(N_Center,1),1)*sDiris.data(j,:)-N_Center;
    Disim_te(j) = sum(diag(sqrt(sub*sub')));
end
Con_lim = control_lim*ones(1,length(Disim_te));
Disim_te = [zeros(1,wid-1),Disim_te];
Con_lim = [zeros(1,wid-1),Con_lim];
fdr = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; 
fprintf('fdr = %.2f\n',fdr)
fdr = roundn(fdr,-2);
far = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  
fprintf('far = %.2f\n',far)
far = roundn(far,-2);
figure(); hold on 
plot([1:length(Disim_te)],Disim_te,'b');
plot([1:length(Disim_te)],control_lim,'r');
title(strcat('FAR:',num2str(far),'%','  FDR:',num2str(fdr),'%'),'FontName','Times','FontSize',13.5);
plot([160,160],[0,max(Disim_te)+10],'k--')


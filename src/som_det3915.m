clear
clc
close all
dnorm = load('D:\Paper 2\te_process\d00.dat');    
d01 = load('D:\Paper 2\te_process\d03.dat');   
dtrain = [dnorm;d01]; 
wid = 20;
[dtrain,Dx] = mncn(dtrain);
for i = 1:size(dtrain,1)-wid+1  D(i,:) = moving_window(dtrain,i+wid-1,wid); end

label_0 = 0*ones(size(dnorm,1)-wid+1,1); label_1 = 1*ones(size(d01,1),1);
label = [label_0;label_1];

sDiris = som_data_struct(D,'name','Te (train)');
sDiris = som_label(sDiris,'add',(find(label==0)),'N');
sDiris = som_label(sDiris,'add',(find(label==1)),'F');
sMap = som_make(sDiris,'small','tracking',1);  
sMap = som_autolabel(sMap,sDiris,'vote');
h0 = som_hits(sMap,sDiris.data(label==0,:));
h1 = som_hits(sMap,sDiris.data(label==1,:));
bmus = som_bmus(sMap,sDiris,1);
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels')
som_show_add('hit',[h0, h1],'MarkerColor',[1 0 0; 0 1 0],'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
for i = 1:size(label,1)
    ind =  sMap.labels(bmus(i));
switch ind{1}
    case 'N',  pre_label(i)  = 0;
        case 'F',  pre_label(i)  = 1;
end
end
pre_label = pre_label';
fdr = sum(label==pre_label)/size(label,1);
title(strcat('Train Result & fdr = ',num2str(fdr)));
figure
plot([1:size(pre_label,1)],pre_label,'ro'); title('Train')
axis([0 1000 -1 2]);

[~,ind] = sort(h0,'descend');
% [~,ind] = sort(h1,'descend');
N_Center = sMap.codebook(ind(1),:);
% 
sub = sMap.codebook(bmus,:)-ones(size(sDiris.data,1),1)*N_Center;
Disim = diag(sqrt(sub*sub'));
[fun,varx] = ksdensity(Disim,'bandwidth',0.2,'function','cdf');
Seq = find(h0>0);
for i = 1:length(Seq)
    D = sDiris.data(bmus==Seq(i),:);  
      d{i} = sMap.codebook(Seq(i),:);
        Vx{i} = D\(ones(size(D,1),1)*d{i});
end
% figure()
% plot(varx,fun,'r-'); title('CDF of Disim');
% control_lim = varx(min(find(fun>0.98)));
clear sDiris
% sMap.codebook(h0<=0,:) = []; 

% 
d01 = load('D:\Paper 2\te_process\d03_te.dat'); 
d01 = scale(d01,Dx(1,:),Dx(2,:));
for i = 1:size(d01,1)-wid+1  Dte(i,:) = moving_window(d01,i+wid-1,wid); end

sDiris = som_data_struct(Dte,'name','Te (test)');
sDiris = som_label(sDiris,'add',(1:151)','N');
sDiris = som_label(sDiris,'add',(152:size(Dte,1))','F');
bmus = som_bmus(sMap,sDiris);
h0 = som_hits(sMap,sDiris.data(1:151,:));
h1 = som_hits(sMap,sDiris.data(152:size(Dte,1),:));
sMap = som_label(sMap, 'clear', [1:size(sMap.labels,1)]');
sMap = som_autolabel(sMap,sDiris,'vote');
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels')
som_show_add('hit',[h0, h1],'MarkerColor',[1 0 0; 0 1 0],'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
for i = 1:size(label,1)
    ind =  sMap.labels(bmus(i));
switch ind{1}
    case 'N',  pre_label(i)  = 0;
        case 'F',  pre_label(i)  = 1;
end
end
pre_label = pre_label';
fdr = sum(label==pre_label)/size(label,1);
title(strcat('Train Result & fdr = ',num2str(fdr)));
figure
plot([1:size(pre_label,1)],pre_label,'ro'); title('Train')
axis([0 1000 -1 2]);


% for i = 1:size(Dte,1)
%     sub = sDiris.data(i,:)*Vx{bmus(i)}-N_Center;
% %    sub = ones(size(N_Center,1),1)*sDiris.data(i,:)-N_Center;
%     Disim_te(i) = sum(diag(sqrt(sub*sub')));
% end
% Disim_te = [zeros(1,wid-1),Disim_te];
% control_lim = [zeros(1,wid-1),control_lim];
% figure(); hold on
% plot([1:length(Disim_te)],Disim_te,'b');
% plot([1:length(Disim_te)],control_lim,'r.-');
% fdr = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; fdr = roundn(fdr,-2);
% far = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  far = roundn(far,-2);
% title(strcat('FAR:',num2str(far),'%','  FDR:',num2str(fdr),'%'),'FontName','Times','FontSize',13.5);

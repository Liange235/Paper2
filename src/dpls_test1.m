clear
close all
clc

%ind = [1:4 7:11 18 23:32 34:41 42:43]; 
ind = [1:52];n = length(ind);
% ind = [1:4];
%********training*****************%
dnorm = load('D:\Paper 2\te_process\d00.dat'); 
d01 = load('D:\Paper 2\te_process\d08.dat');  d01 = d01(:,ind);
d02 = load('D:\Paper 2\te_process\d09.dat');  d02 = d02(:,ind);
d03 = load('D:\Paper 2\te_process\d10.dat');  d03 = d03(:,ind);
[N_Center,Dx] = mysom_train(1,0.9);
d01 = scale(d01,Dx(1,:),Dx(2,:)); d02 = scale(d02,Dx(1,:),Dx(2,:)); d03 = scale(d03,Dx(1,:),Dx(2,:));

d01 = d01-N_Center(ones(size(d01,1),1),:);d02 = d02-N_Center(ones(size(d02,1),1),:);d03 = d03-N_Center(ones(size(d03,1),1),:);
% [D,Dx] = mncn(D);
d01 = d01+sqrt(diag(d01*d01'))*ones(1,n); d02 = d02+sqrt(diag(d02*d02'))*ones(1,n); d03 = d03+sqrt(diag(d03*d03'))*ones(1,n);
dlim = 480; tao = 51; dr = 0;
d01 = d01(-dr+1:dlim-dr,:);   d02 = d02(-dr+1:dlim-dr,:);   d03 = d03(-dr+1:dlim-dr,:);
% label_1 = 1*ones(size(d01,1),1);  label_2 = 2*ones(size(d02,1),1);  label_3 = 3*ones(size(d03,1),1); %%1 careful
d01 = dynamic(d01,tao);  d02 = dynamic(d02,tao);  d03 = dynamic(d03,tao);
label_1 = 1*ones(size(d01,1),1);  label_2 = 2*ones(size(d02,1),1);  label_3 = 3*ones(size(d03,1),1);
label = [label_1;label_2;label_3];
D_y1 = [label_1 zeros(size(label_1,1),2)]; D_y2 = [zeros(size(label_2,1),1) 1/2*label_2 zeros(size(label_2,1),1)]; D_y3 = [zeros(size(label_3,1),2) 1/3*label_3];
% D_y1 = dynamic(D_y1,tao); D_y2 = dynamic(D_y2,tao); D_y3 = dynamic(D_y3,tao);
% label_1 = 1*ones(size(d01,1),1);  label_2 = 2*ones(size(d02,1),1);  label_3 = 3*ones(size(d03,1),1);
D = [d01;d02;d03];
D_y = [D_y1;D_y2;D_y3];
% [D_y,Dy] = mncn(D_y,1);
% [plsss,cplsss,mlv,bpls,~,loading_x,loading_y] = plscvblk(D,D_y,10,size(d02,2),1);
% mlv = size(d01,2);
mlv = 4;
[p,q,w,t,u,bb,ssqdif] = pls(D,D_y,mlv); 
[bb,ind] = sort(bb,'descend');
% figure
% plot(bb,'k*-');
w = w(:,ind);  p = p(:,ind); q = q(:,ind);
for i = 1:size(p,2)
    t(:,i) = D*w(:,i);
    D = D-t(:,i)*p(:,i)';
end
t = t.*bb(ones(size(D,1),1),:); 
% coef = conpred1(bb,w,p,q,mlv);
% t = D*coef;
% t = rescale(t,Dy);

figure
hold on
plot(t(1:size(label_1,1),1),t(1:size(label_1,1),2),'ro');
plot(t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),1),t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),2),'go');
plot(t(size(label_1,1)+size(label_2,1)+1:end,1),t(size(label_1,1)+size(label_2,1)+1:end,2),'bo');
%********使用SOM网络测试其特征提取效果*****************%
for i = 1:size(t,2) comp{i} = num2str(i); end
 sDiris = som_data_struct([t],'name','Te (train)',...
			  'comp_names',comp);
sDiris = som_label(sDiris,'add',(1:size(label_1,1))','f1');
sDiris = som_label(sDiris,'add',(size(label_1,1)+1:size(label_1,1)+size(label_2,1))','f2');
sDiris = som_label(sDiris,'add',(size(label_1,1)+size(label_2,1)+1:size(D,1))','f3');
sMap = som_make(sDiris,'small','tracking',1); 
sMap = som_autolabel(sMap,sDiris,'vote');
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels')
h1 = som_hits(sMap,sDiris.data(1:size(label_1,1),:));
h2 = som_hits(sMap,sDiris.data(size(label_1,1)+1:size(label_1,1)+size(label_2,1),:));
h3 = som_hits(sMap,sDiris.data(size(label_1,1)+size(label_2,1)+1:size(t,1),:));
som_show_add('hit',[h1, h2, h3],'MarkerColor',[1 0 0; 0 1 0; 0 0 1],'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
bmus = som_bmus(sMap,sDiris,1);
for i = 1:size(label,1)
    ind =  sMap.labels(bmus(i));
switch ind{1}
    case 'f1',  pre_label(i)  = 1;
        case 'f2',  pre_label(i)  = 2;
            case 'f3',  pre_label(i)  = 3;
end
end
pre_label = pre_label';
fdr = sum(label==pre_label)/size(label,1);
title(strcat('Train Result & fdr = ',num2str(fdr)));
figure
plot([1:size(pre_label,1)],pre_label,'ro'); title('Train')
clear sDiris pre_label t

%%%***************利用som得到每类数据的聚类中心****************%%%%%%
% for i = 1:3
% sDiris_dum = som_data_struct(D(label==i,:));
% sMap_dum = som_make(sDiris_dum,'tracking',0);
% [bmus,~]= som_bmus(sMap_dum,D(label==i,:),1); %算出每个观测向量D对应的1-levelBMU
% for j = 1:max(bmus) freq(j) = sum(bmus==j); end  
% [~,ind] = sort(freq);
% Center(i,:) = sMap_dum.codebook(ind(end),:);  % 样本的聚类中心
% clear freq
% end
%%****************************************************************************************%
%%********testing*****************%
ind = [1:52];
% ind = [1:4];
d01_te = load('D:\Paper 2\te_process\d08_te.dat');  d01_te = d01_te(:,ind);
d02_te = load('D:\Paper 2\te_process\d09_te.dat');  d02_te = d02_te(:,ind);
d03_te = load('D:\Paper 2\te_process\d10_te.dat');  d03_te = d03_te(:,ind);
d01_te = scale(d01_te,Dx(1,:),Dx(2,:)); d02_te = scale(d02_te,Dx(1,:),Dx(2,:)); d03_te = scale(d03_te,Dx(1,:),Dx(2,:));
dlim = 200;  dr = 0;
% Dte = scale(Dte,Dx(1,:),Dx(2,:));
d01_te = d01_te-N_Center(ones(size(d01_te,1),1),:);d02_te = d02_te-N_Center(ones(size(d02_te,1),1),:);d03_te = d03_te-N_Center(ones(size(d03_te,1),1),:);
% d01_te = Dte(1:size(d01_te,1),:); d02_te = Dte(size(d01_te,1)+1:size(d01_te,1)+size(d02_te,1),:); d03_te = Dte(size(d01_te,1)+size(d02_te,1)+1:end,:);
d01_te = d01_te+sqrt(diag(d01_te*d01_te'))*ones(1,n); d02_te = d02_te+sqrt(diag(d02_te*d02_te'))*ones(1,n); d03_te = d03_te+sqrt(diag(d03_te*d03_te'))*ones(1,n);
d01_te = d01_te(160-dr+1:160+dlim-dr,:);   d02_te = d02_te(160-dr+1:160+dlim-dr,:);   d03_te = d03_te(160-dr+1:160+dlim-dr,:);

%* Dte = scale(Dte,Dx(1,:),Dx(2,:));  d01_te = Dte(1:size(d01_te,1),:); d02_te = Dte(size(d01_te,1)+1:size(d01_te,1)+size(d02_te,1),:); d03_te = Dte(size(d01_te,1)+size(d02_te,1)+1:end,:);
d01_te = dynamic(d01_te,tao);  d02_te = dynamic(d02_te,tao);  d03_te = dynamic(d03_te,tao);
Dte = [d01_te;d02_te;d03_te];
label_1 = 1*ones(size(d01_te,1),1);  label_2 = 2*ones(size(d02_te,1),1);  label_3 = 3*ones(size(d03_te,1),1);
label = [label_1;label_2;label_3];
%%%%%%此处只考虑两种情况的数据：正常+异常态数据
%%%计算测试数据对于正常和
% sum_de1 = sum((Dte(label==1,:)-ones(size(d01_te,1),1)*Center(1,:)).^2,2);
% sum_de2 = sum((Dte(label==2,:)-ones(size(d02_te,1),1)*Center(1,:)).^2,2);
% sum_de3 = sum((Dte(label==3,:)-ones(size(d03_te,1),1)*Center(1,:)).^2,2);
% mu1 = (sum_de1)./(sum_de1+sum_de2+sum_de3); mu1 = mu1.^-1; 
% mu2 = (sum_de2)./(sum_de1+sum_de2+sum_de3); mu2 = mu2.^-1;
% mu3 = (sum_de3)./(sum_de1+sum_de2+sum_de3); mu3 = mu3.^-1;
% figure()
% hold on 
% plot(mu3,'r'); %plot(mu2,'g');  plot(mu3,'b');  
% plot([161,161],[min(mu1),max(mu1)],'k-');
%legend('normal center','anomalous center');
% legend('fault1','fault2','fault3');
% [f1,xi1,w] = ksdensity(mu1(1:160,:),'support','positive','npoints',1000);
% [f2,xi2,w] = ksdensity(mu1(161:480,:),'support','positive','npoints',1000);
% [f3,xi3,w] = ksdensity(mu3,'support','positive','npoints',1000);
% figure()
% hold on
% plot(xi1,f1,'r'); plot(xi2,f2,'g'); plot(xi3,f3,'b')
% legend('normal center','anomalous center');
% legend('fault1','fault2','fault3'); title('CDF of three Time-Series Data')
%%%**************************
% t = Dte*w.*bb(ones(size(Dte,1),1),:);  %%%%%%关键一句话
for i = 1:size(p,2)
    t(:,i) = Dte*w(:,i);
    Dte = Dte-t(:,i)*p(:,i)';
end
t = t.*(bb(ones(size(Dte,1),1),:));
% t = Dte*coef;
% t = rescale(t,Dy);
figure(); plot(t); legend('1','2','3');
figure
hold on
plot(t(1:size(label_1,1),1),t(1:size(label_1,1),2),'ro');
plot(t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),1),t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),2),'go');
plot(t(size(label_1,1)+size(label_2,1)+1:end,1),t(size(label_1,1)+size(label_2,1)+1:end,2),'bo');
for i = 1:size(t,2) comp{i} = num2str(i); end
 sDiris = som_data_struct([t],'name','Te (test)',...
			  'comp_names',comp);
sDiris = som_label(sDiris,'add',(1:size(label_1,1))','f1');
sDiris = som_label(sDiris,'add',(size(label_1,1)+1:size(label_1,1)+size(label_2,1))','f2');
sDiris = som_label(sDiris,'add',(size(label_1,1)+size(label_2,1)+1:size(t,1))','f3');
% sMap = som_make(sDiris,'small','tracking',1); 
% sMap = som_label(sMap, 'clear', [1:size(sMap.labels,1)]');
% sMap = som_autolabel(sMap,sDiris,'vote');
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels');


% sMap.codebook((h1+h2+h3)<1,:) = [];
% sMap.labels((h1+h2+h3)<1,:) = [];

h1_te = som_hits(sMap,sDiris.data(1:size(label_1,1),:));
h2_te = som_hits(sMap,sDiris.data(size(label_1,1)+1:size(label_1,1)+size(label_2,1),:));
h3_te = som_hits(sMap,sDiris.data(size(label_1,1)+size(label_2,1)+1:size(t,1),:));
som_show_add('hit',[h1_te, h2_te, h3_te],'MarkerColor',[1 0 0; 0 1 0; 0 0 1],'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
pre_label = zeros(size(label));
bmus = som_bmus(sMap,sDiris,1);
for i = 1:size(label,1)
    ind =  sMap.labels(bmus(i));
switch ind{1}
    case 'f1',  pre_label(i)  = 1;
        case 'f2',  pre_label(i)  = 2;
            case 'f3',  pre_label(i)  = 3;
end
end
fdr = sum(label==pre_label)/size(label,1);
title(strcat('Test Result & fdr = ',num2str(fdr)));
% [mqe,tge,cbe] = som_quality(sMap,sDiris); cbe
figure; hold on
plot([1:size(pre_label,1)],pre_label,'ro'); title('Test')
plot([max(find(label==1)),max(find(label==1))],[1,2],'k--','LineWidth',1);plot([max(find(label==2)),max(find(label==2))],[2,3],'k--','LineWidth',1);
% lg = [h1,h1_te,h2,h2_te,h3,h3_te];
% figure()
% bar([h3,h3_te]);
% set(gca,'xtick',[1:length(h1)]);
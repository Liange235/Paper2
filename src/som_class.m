clear 
clc
close all
% global shift
%% data preparation
D = []; label = [];hits = []; lg = 0;
lib = [0,8,10]'; numel = length(lib);
dlim = 480*ones(1,numel);
tao = 30; dr = 00*ones(1,numel); 
[N_Center,Dx] = mysom_train(1,0.9); 
Sf = cell(1,numel+1);
Sf(1) = {''};
for i = 2:numel+1
     str = strcat('f',num2str(lib(i-1)));
     Sf(i) = cellstr(str);
end
for i = 1:numel
str = num2str(lib(i));
if (lib(i)<10)
    str = strcat('0',num2str(lib(i)));
end
d_com = load(strcat('D:\研究生\Paper 2\te_process\d',str,'.dat')) ;
d_com = d_com(1:480,:);
d_com = scale(d_com,Dx(1,:),Dx(2,:));
absd_com = zeros(480,52);
% absd_com = [];
for j = 1:lg
    sub = N_Center(j,:);
    stat = sqrt(diag((d_com-sub(ones(size(d_com,1),1),:))*(d_com-sub(ones(size(d_com,1),1),:))'));
    absd_com = absd_com+stat*ones(1,52);
%     absd_com = [absd_com,stat];
    d_com = d_com-ones(size(d_com,1),1)*sub;
end
d_com = d_com+absd_com;
d_com = d_com(dr(i)+1:dlim(i)+dr(i),:);
dd_com = dynamic(d_com,tao);
D = [D;dd_com];
label_com = i*ones(size(dd_com,1),1);
label = [label;label_com];
end
D_y = zeros(size(D,1),numel); 
for i = 1:numel D_y(label==i,i)=1*ones(sum(label==i),1);  end                                                                                     
% Wxpca = mypca(D);
% D = D*Wxpca;
% [D,Sd] = mncn(D);
% coef = canonize(D,label);
% t = D*coef(:,1:numel-1);
%% pls feature extraction
% [plsss,cplsss,mlv,bpls,~,loading_x,loading_y] = plscvblk(D,D_y,10,size(dd_com,2),1);
%%****pls特征提取******%%%%
mlv = 25;
[p,q,w,t,u,bb,ssqdif] = pls(D,D_y,mlv);
[bb,ind] = sort(bb,'descend');
w = w(:,ind);  p = p(:,ind); q = q(:,ind); 
for i = 1:size(w,2)
    t(:,i) = D*w(:,i);
    D = D-t(:,i)*p(:,i)';
end
t = t.*bb(ones(size(D,1),1),:);
%%%***********************
% coef = conpred1(bb,w,p,q,mlv);
% t = D*coef;
% t = rescale(t,Dy);
%% som clustering
sDiris = som_data_struct(t,'name','Te (train)');
for i = 1:numel
     sDiris = som_label(sDiris,'add',find(label==i),Sf(i+1));
end
sMap = som_make(sDiris,'small','tracking',1); 
sMap = som_autolabel(sMap,sDiris,'vote');
bmus = som_bmus(sMap,sDiris,1);
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels')
for i = 1:numel
    h = som_hits(sMap,sDiris.data(label==i,:));
    hits = [hits,h];
end
for i = 1:size(label,1)
    ind =  sMap.labels(bmus(i));
   pre_label(i) = find(strcmp(ind{1},Sf))-1;
end
pre_label = pre_label';
% calculate the overall classification rate (training)
fdr = sum(label==pre_label)/size(label,1);
% calculate the classification rate of each fault (training)
for i = 1:numel
    label_com = label(label==i);
    plabel_com = pre_label(label==i);
%     fdr_com(i) = sum(label_com==plabel_com)/size(label_com,1);
    tp = sum(label_com==plabel_com); fp = sum(pre_label(label~=i)==i);
    precision = tp/(tp+fp);
    recall = tp/size(label_com,1);
    fdr_com(i) = 2*(precision*recall)/(precision+recall);
end
% Fs_m = mean(fdr_com);
color = hsv;
som_show_add('hit',hits,'MarkerColor',color(round(linspace(1,54,numel)),:),'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
title(strcat('Offline Result & fdr = ',num2str(fdr)));
figure; hold on
plot(pre_label,'r.');
subs = dlim-tao; cumx(1) = subs(1);
plot([1,cumx(1)],[1,1],'k--','LineWidth',2);%%*
plot([cumx(1),cumx(1)],[0,numel],'k--','LineWidth',1);
text(1/2*cumx(1),1+1/3,...
	['F-s = ',num2str(fdr_com(1))],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],'FontName','Times');
for i = 2:numel
    cumx(i) = sum(subs(1:i));
    plot([cumx(i-1)+1,cumx(i)],[i,i],'k--','LineWidth',2); %%*
    plot([cumx(i),cumx(i)],[0,numel],'k--','LineWidth',1);
    text(cumx(i-1)+1/2*(cumx(i)-cumx(i-1)),1/3+i,...
	['F-s = ',num2str(fdr_com(i))],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],'FontName','Times');
end
axis([0,size(label,1),0,numel+1/2]);
set(gca,'YTick',[0:numel]); set(gca,'YTickLabel',Sf);
set(gca,'xTick',[0,cumx]);
title('Offline development result','FontName','Times','FontSize',14)
xlabel('Samples','FontName','Times','Fontsize',14);
ylabel('Type','FontName','Times','Fontsize',14);

clear sDiris t pre_label dlim subs label_com
% sMap.codebook(sum(hits,2)<fil_n,:) = [];

%%%%%%***********test**************%%%%%%%%%%%%%
Dte = []; label_te = []; hits_te = []; 
dlim = 200*ones(1,numel);  
% dr =  [0,0,0,50,0,0,52,64];
% dr = [50,50,50,72];
for i = 1:numel
str = num2str(lib(i));  
if (lib(i)<10)
    str = strcat('0',num2str(lib(i)));
end
d_com = load(strcat('D:\研究生\Paper 2\te_process\d',str,'_te.dat'));
d_com = scale(d_com,Dx(1,:),Dx(2,:));
absd_com = zeros(960,52);
% absd_com = [];
for j = 1:lg
    sub = N_Center(j,:);
    stat = sqrt(diag((d_com-sub(ones(size(d_com,1),1),:))*(d_com-sub(ones(size(d_com,1),1),:))'));
    absd_com = absd_com+stat*ones(1,52);
%      absd_com = [absd_com,stat];
    d_com = d_com-ones(size(d_com,1),1)*sub;
end
d_com = d_com+absd_com; 
d_com = d_com(160+dr(i)+1:160+dlim(i)+dr(i),:);
dd_com = dynamic(d_com,tao);
Dte = [Dte;dd_com];
label_com = i*ones(size(dd_com,1),1);
label_te = [label_te;label_com];
end
% Dte = scale(Dte,Sd(1,:),Sd(2,:));
% t = Dte*coef(:,1:numel-1);
% Dte = Dte*Wxpca;
%%****pls特征提取******%%%%
for i = 1:size(w,2)
    t(:,i) = Dte*w(:,i);
    Dte = Dte-t(:,i)*p(:,i)';
end
t = t.*bb(ones(size(Dte,1),1),:);%% pls projection matrix*Dtest
%%%***********************
% t = Dte*coef;
% t = rescale(t,Dy);
sDiris = som_data_struct(t,'name','Te (test)');
% sMap = som_make(sDiris,'small','tracking',1);
% sMap = som_label(sMap, 'clear', [1:size(sMap.labels,1)]');
% sMap = som_autolabel(sMap,sDiris,'vote');
%*******som bmu filtering**********%
for i = 1:numel
    h = som_hits(sMap,sDiris.data(label_te==i,:));
    hits_te = [hits_te,h];  
end
% sMap.codebook(sum(hits_te,2)<fil_n,:) = [];
% sMap.labels(sum(hits_te,2)<fil_n,:) = [];
bmus = som_bmus(sMap,sDiris,1);

for i = 1:size(label_te,1)
    ind =  sMap.labels(bmus(i));
    pre_label(i) = find(strcmp(ind{1},Sf))-1;
end
pre_label = pre_label';

% calculate the overall classification rate (testing)
fdr = sum(label_te==pre_label)/size(label_te,1);
% calculate the classification rate of each fault (testing)
for i = 1:numel
    label_com = label_te(label_te==i);
    plabel_com = pre_label(label_te==i);
%     fdr_com(i) = sum(label_com==plabel_com)/size(label_com,1);
    tp = sum(label_com==plabel_com); fp = sum(pre_label(label_te~=i)==i);
    precision = tp/(tp+fp);
    recall = tp/size(label_com,1);
    fdr_com(i) = 2*(precision*recall)/(precision+recall);
end
% Fs_m = mean(fdr_com);
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels');
som_show_add('hit',hits_te,'MarkerColor',color(round(linspace(1,54,numel)),:),'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','b','Subplot',2);
title(strcat('Online Result & fdr = ',num2str(fdr)));
figure; hold on
plot(pre_label,'r.');
subs = dlim-tao; cumx(1) = subs(1);
plot([1,cumx(1)],[1,1],'k--','LineWidth',2); %%*
plot([cumx(1),cumx(1)],[0,numel],'k--','LineWidth',1);
text(1/2*cumx(1),1+1/3,...
	['F-s = ',num2str(fdr_com(1))],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],'FontName','Times');
for i = 2:numel
    cumx(i) = sum(subs(1:i));
    plot([cumx(i-1)+1,cumx(i)],[i,i],'k--','LineWidth',2); %%*
    plot([cumx(i),cumx(i)],[0,numel],'k--','LineWidth',1);
    text(cumx(i-1)+1/2*(cumx(i)-cumx(i-1)),1/3+i,...
	['F-s = ',num2str(fdr_com(i))],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],'FontName','Times');
end
axis([0,size(label_te,1),0,numel+1/2]);
set(gca,'YTick',[0:numel]); set(gca,'YTickLabel',Sf);
set(gca,'xTick',[0,cumx]);
title('Online implementation result','FontName','Times','FontSize',14)
xlabel('Samples','FontName','Times','Fontsize',14);
ylabel('Type','FontName','Times','Fontsize',14);

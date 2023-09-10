close all
clear all
clc
%--------------------------------------------------
D = []; label = [];lg = 0;tao = 0;
lib = [0,1,2,4,5,6,7,8,10,11,12,14,17,18,20]'; numel = length(lib);
% [N_Center,Dx] = mysom_train(1,0.9); 
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
d_com = load(strcat('E:\My_Education\Graduate\Paper 2\te_process\d',str,'.dat')) ;
d_com = d_com(1:480,:);
% d_com = scale(d_com,Dx(1,:),Dx(2,:));
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
d_com = d_com(1:240,:);
dd_com = dynamic(d_com,tao);
D = [D;dd_com];
label_com = i*ones(size(dd_com,1),1);
label = [label;label_com];
end
[D,Sd] = mncn(D);

%-----------------------------------------------------
%%   SVM 模型Learning and Learning Parameters
c = 1000;
lambda = 1e-7;
kerneloption= 2;
kernel='gaussian';
verbose = 1;
[xsup,w,b,nbsv]=svmmulticlassoneagainstall(D,label,numel,c,lambda,kernel,kerneloption,verbose);
%%
%%%%%%***********test**************%%%%%%%%%%%%%
Dte = []; label_te = []; 
for i = 1:numel
str = num2str(lib(i));  
if (lib(i)<10)
    str = strcat('0',num2str(lib(i)));
end
d_com = load(strcat('E:\My_Education\Graduate\Paper 2\te_process\d',str,'_te.dat'));
% d_com = scale(d_com,Dx(1,:),Dx(2,:));
absd_com = zeros(960,52);
for j = 1:lg
    sub = N_Center(j,:);
    stat = sqrt(diag((d_com-sub(ones(size(d_com,1),1),:))*(d_com-sub(ones(size(d_com,1),1),:))'));
    absd_com = absd_com+stat*ones(1,52);
%      absd_com = [absd_com,stat];
    d_com = d_com-ones(size(d_com,1),1)*sub;
end
d_com = d_com+absd_com; 
d_com = d_com(160+1:end,:);
dd_com = dynamic(d_com,tao);
Dte = [Dte;dd_com];
label_com = i*ones(size(dd_com,1),1);
label_te = [label_te;label_com];
end
Dte = scale(Dte,Sd(1,:),Sd(2,:));
[pre_label,maxi] = svmmultival(Dte,xsup,w,b,nbsv,kernel,kerneloption);
for i = 1:numel
    label_com = label_te(label_te==i);
    plabel_com = pre_label(label_te==i);
%     fdr_com(i) = sum(label_com==plabel_com)/size(label_com,1);
    tp = sum(label_com==plabel_com); fp = sum(pre_label(label_te~=i)==i);
    precision = tp/(tp+fp);
    recall = tp/size(label_com,1);
    fdr_com(i) = 2*(precision*recall)/(precision+recall);
end
figure; hold on
plot(pre_label,'r.');
subs = 800*ones(1,numel); ; cumx(1) = subs(1);
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
print(gcf,'-dbitmap','E:\My_Education\Graduate\Paper 2\论文\Revision\Pic\svm2.bmp');
mail2me('title','done','E:\My_Education\Graduate\Paper 2\论文\Revision\Pic\svm2.bmp');

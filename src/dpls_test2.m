clear
close all
clc
%********training*****************%
d01 = load('D:\Paper 2\te_process\d05.dat');  
d02 = load('D:\Paper 2\te_process\d04.dat');  
d03 = load('D:\Paper 2\te_process\d02.dat');  
% [N_Center,Dx] = mysom_train(1,0.9);
% d01 = scale(d01,Dx(1,:),Dx(2,:)); d02 = scale(d02,Dx(1,:),Dx(2,:)); d03 = scale(d03,Dx(1,:),Dx(2,:));
dlim = 480; tao = 0; dr = 0;
d01 = d01(-dr+1:dlim-dr,:);   d02 = d02(-dr+1:dlim-dr,:);   d03 = d03(-dr+1:dlim-dr,:);
d01 = dynamic(d01,tao);  d02 = dynamic(d02,tao);  d03 = dynamic(d03,tao);
label_1 = 1*ones(size(d01,1),1);  label_2 = 2*ones(size(d02,1),1);  label_3 = 3*ones(size(d03,1),1);
label = [label_1;label_2;label_3];
D_y1 = [label_1 zeros(size(label_1,1),2)]; D_y2 = [zeros(size(label_2,1),1) 1/2*label_2 zeros(size(label_2,1),1)]; D_y3 = [zeros(size(label_3,1),2) 1/3*label_3];
D = [d01;d02;d03];
[D,Dx] = mncn(D);
D_y = [D_y1;D_y2;D_y3];
% [plsss,cplsss,mlv,bpls,~,loading_x,loading_y] = plscvblk(D,D_y,10,size(d02,2),1);
mlv = 2;
[p,q,w,t,u,bb,ssqdif,p1] = pls(D,D_y,mlv);
coef = conpred1(bb,w,p,q,mlv);
t = D*coef;
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
clear sDiris pre_label

%%%%%**************test*************%%%%%%%%%
d01_te = load('D:\Paper 2\te_process\d05_te.dat'); 
d02_te = load('D:\Paper 2\te_process\d04_te.dat');  
d03_te = load('D:\Paper 2\te_process\d02_te.dat');  
% d01_te = scale(d01_te,Dx(1,:),Dx(2,:)); d02_te = scale(d02_te,Dx(1,:),Dx(2,:)); d03_te = scale(d03_te,Dx(1,:),Dx(2,:));
d01_te = d01_te(160-dr+1:160+dlim-dr,:);   d02_te = d02_te(160-dr+1:160+dlim-dr,:);   d03_te = d03_te(160-dr+1:160+dlim-dr,:);
d01_te = dynamic(d01_te,tao);  d02_te = dynamic(d02_te,tao);  d03_te = dynamic(d03_te,tao);
Dte = [d01_te;d02_te;d03_te];
Dte = scale(Dte,Dx(1,:),Dx(2,:));
label_1 = 1*ones(size(d01_te,1),1);  label_2 = 2*ones(size(d02_te,1),1);  label_3 = 3*ones(size(d03_te,1),1);
label = [label_1;label_2;label_3];
t = Dte*coef;
for i = 1:size(t,2) comp{i} = num2str(i); end
 sDiris = som_data_struct([t],'name','Te (test)',...
			  'comp_names',comp);
sDiris = som_label(sDiris,'add',(1:size(label_1,1))','f1');
sDiris = som_label(sDiris,'add',(size(label_1,1)+1:size(label_1,1)+size(label_2,1))','f2');
sDiris = som_label(sDiris,'add',(size(label_1,1)+size(label_2,1)+1:size(t,1))','f3');
figure
colormap(1-gray)
som_show(sMap,'umat',{'all','U-matrix'},'empty','Labels');
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
figure
plot([1:size(pre_label,1)],pre_label,'ro'); title('Test')
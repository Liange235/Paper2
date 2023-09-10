clear
close all
clc
%********training*****************%
d00 = load('D:\Paper 2\te_process\d00.dat');
d01 = load('D:\Paper 2\te_process\d04.dat');
d02 = load('D:\Paper 2\te_process\d05.dat');
d = d00(1:20,:);
t = 15;
% d01 = [d;d01];
% d02 = [d;d02];
[d00 d00x] = mncn(d00);  [d01 d01x] = mncn(d01);  [d02 d02x] = mncn(d02);
% d00 = dynamic(d00,t);  d01 = dynamic(d01,t);  d02 = dynamic(d02,t);
D = [d00;d01;d02];
label_1 = 1*ones(size(d00,1),1);  label_2 = 100*ones(size(d01,1),1);  label_3 = 210*ones(size(d02,1),1);
label = [label_1;label_2;label_3];
D_y = [label_1 zeros(size(label_1,1),2);zeros(size(label_2,1),1) 1/2*label_2 zeros(size(label_2,1),1);zeros(size(label_3,1),2) 1/3*label_3];
D_y = mncn(D_y);
% [plsss,cplsss,mlv,bpls,~,loading_x,loading_y] = plscvblk(D,D_y,8,10,1);
figure(1)
hold on
[p,q,w,t,u,bb,ssqdif] = pls(D,D_y,40);
b = conpred1(bb,w,p,q,40);
% t = D*b;
plot(t(1:size(label_1,1),1),t(1:size(label_1,1),2),'bo');
plot(t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),1),t(size(label_1,1)+1:size(label_1,1)+size(label_2,1),2),'ro');
plot(t(size(label_1,1)+size(label_2,1)+1:end,1),t(size(label_1,1)+size(label_2,1)+1:end,2),'go');
%********使用SOM网络测试其特征提取效果*****************%
for i = 1:size(t,2) comp{i} = num2str(i); end
 sDiris = som_data_struct([t],'name','Te (simulated)',...
			  'comp_names',comp);
sDiris = som_label(sDiris,'add',(1:size(label_1,1))','cluster1');
sDiris = som_label(sDiris,'add',(size(label_1,1)+1:size(label_1,1)+size(label_2,1))','cluster2');
sDiris = som_label(sDiris,'add',(size(label_1,1)+size(label_2,1)+1:size(t,1))','cluster3');
sMap = som_make(sDiris);
sMap = som_autolabel(sMap,sDiris,'vote');
figure(2)
colormap(1-gray)
som_show(sMap,'umat','all','empty','Labels')
som_show_clear('hit',1)
h1 = som_hits(sMap,sDiris.data(1:size(label_1,1),:));
h2 = som_hits(sMap,sDiris.data(size(label_1,1)+1:size(label_1,1)+size(label_2,1),:));
h3 = som_hits(sMap,sDiris.data(size(label_1,1)+size(label_2,1)+1:size(t,1),:));
som_show_add('hit',[h1, h2, h3],'MarkerColor',[1 0 0; 0 1 0; 0 0 1],'Subplot',1);
som_show_add('label',sMap,'Textsize',8,'TextColor','r','Subplot',2);
% %********testing*****************%
% d00_te = load('D:\研究生\Paper2\te_process\d00_te.dat');
% d01_te = load('D:\研究生\Paper2\te_process\d01_te.dat');
% d02_te = load('D:\研究生\Paper2\te_process\d02_te.dat');


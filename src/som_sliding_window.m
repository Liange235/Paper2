clear
clc
close all
par = [0.95,0.99,0.994,0.998,0.85,0.95,0.99,0.999,0.999,0.9,0.96,0.99,0.99,0.99,0.99,0.999,0.93,0.97,0.97,0.84,0.999];
wid = [1,3,1,2,16,1,1,1,1,8,4,1,2,1,2,1,1,2,1,1,1];
for i = 1:21
[N_Center,Dx,control_lim] = mysom_train(wid(i),par(i));
str = num2str(i);
if (i<10)
    str = strcat('0',num2str(i));
end
d01 = load(strcat('D:\Paper 2\te_process\d',str,'_te.dat'));
for j = 1:size(d01,1)-wid(i)+1  Dte(j,:) = moving_window(d01,j+wid(i)-1,wid(i)); end
Dte = scale(Dte,Dx(1,:),Dx(2,:));
sDiris = som_data_struct(Dte,'name','Te (test)');
% h1 = som_hits(sMap,sDiris.data);
% bmus = som_bmus(sMap,sDiris);
for j = 1:size(Dte,1)
%       sub = sDiris.data(i,:)*Vx{bmus(i)}-norm_center;
%       Disim_te(i) = norm(sub);
%       Disim_te(i) = norm(sDiris.data(i,:)-norm_center);
        sub = ones(size(N_Center,1),1)*sDiris.data(j,:)-N_Center;
        Disim_te(j) = sum(diag(sqrt(sub*sub')));
end
Con_lim = control_lim*ones(1,length(Disim_te));
Disim_te = [zeros(1,wid(i)-1),Disim_te];
Con_lim = [zeros(1,wid(i)-1),Con_lim];
int_var = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; fdr(i) = roundn(int_var,-2);
int_var = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  far(i) = roundn(int_var,-2);
D_str{i} = [Disim_te',Con_lim'];
clear Disim_te Dte
end
multi_plot(D_str);
figure(); hold on 
x = [1:21]; y = [fdr;-far]';
bar(x,y);
sub = legend('Fdr','Far');
set(sub,'FontName','Times','FontSize',13.5);
set(gca,'xtick',[1:22]); plot([x,22],-5*ones(1,22),'k--');plot([x,22],-10*ones(1,22),'k--')
set(gca,'xlim',[0,22]);
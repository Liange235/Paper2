%%%% moving window based on som for fault detection
clear
clc
close all
dnorm = load('D:\研究生\Paper 2\te_process\d00.dat'); D = dnorm;
% d01 = load('D:\Paper 2\te_process\d04.dat'); d01 = d01(:,ind);
%label_1 = 1*ones(size(d01,1),1);
wid = 1;%*% window width

% for i = 1:size(d,1)-wid+1  D(i,:) = moving_window(d,i+wid-1,wid); end
[D,Dx] = mncn(D);
label_0 = 0*ones(size(D,1),1);  label = label_0;
for i = 1:size(D,2) comp{i} = num2str(i); end
 sDiris = som_data_struct([D],'name','Te (train)',...
			  'comp_names',comp);
sDiris = som_label(sDiris,'add',(find(label==0))','norm');
sMap = som_make(sDiris,'small','tracking',1);  
sMap = som_autolabel(sMap,sDiris,'vote');

h0 = som_hits(sMap,sDiris.data(label==0,:));
[~,ind] = sort(h0,'descend');
N_Center = sMap.codebook(ind(1),:); 
bmus = som_bmus(sMap,sDiris.data(label==0,:),1);
% sub = sMap.codebook(bmus,:)-ones(size(bmus,1),1)*norm_center(1,:);
% Disim = diag(sqrt(sub*sub'));
sub = sDiris.data-ones(size(bmus,1),1)*N_Center;
Disim = diag(sqrt(sub*sub'));
% for i = 2:28
%     N_Center(i,:) = sMap.codebook(ind(i),:);
%     sub = sDiris.data-ones(size(bmus,1),1)*N_Center(i,:);
%     Disim = Disim+diag(sqrt(sub*sub'));
% end

[fun,varx] = ksdensity(Disim,'bandwidth',0.1,'function','cdf');
% figure()
% plot(varx,fun,'r-'); title('CDF of Disim');
control_lim = varx(min(find(fun>0.9)));
% Seq = find(h0>5);
% for i = 1:length(Seq)
%     D = sDiris.data(bmus==Seq(i),:);  
%       d = sMap.codebook(Seq(i),:);
%         Vx{i} = D\(d(ones(size(D,1),1),:));
% end

% sMap.codebook(h0<=10,:) = []; 

% figure
% colormap(1-gray)
% som_show(sMap,'umat','all','empty','labels');
% som_show_add('hit',[h0,h1],'MarkerColor',[1 0 0;0 1 0],'markersize',0.5,'subplot',1);
% som_show_add('label',sMap,'Textsize',8,'TextColor','b','subplot',2);
% sMap = smoothing(sMap,h1);
% h0 = som_hits(sMap,sDiris.data(label==0,:));
% h1 = som_hits(sMap,sDiris.data(label==1,:));
% figure
% colormap(1-gray)
% som_show(sMap,'umat','all')
% som_show_add('hit',[h0,h1],'MarkerColor',[1 0 0;0 1 0]);
% bmus = som_bmus(sMap,sDiris,1);
% center_norm = sMap.codebook(ind(1),:);
% U = som_umat(sMap);
% Um = U(1:2:size(U,1),1:2:size(U,2));
% %Um = 1-Um(:)/max(Um(:)); Um(find(h1==0)) = 0;
% figure(); colormap(1-gray)
% % code = som_colorcode(sMap,'rgb1');
% h = som_cplane(sMap,Um(:));
% set(h,'Edgecolor','none');
% hold on
% som_grid(sMap,'Label',cellstr(int2str(h0)),...
% 	 'Line','none','Marker','none','Labelcolor','m');
% set(gcf, 'InvertHardCopy', 'off');
clear sDiris
%***************Test the fault detection rate of SOM method*************%
% ind = [1:4 7:11 18 23:32 34:41 42 43];
par = [0.95,0.95,0.994,0.93,0.9,0.95,0.99,0.999,0.999,0.95,0.95,0.99,0.99,0.99,0.99,0.999,0.93,0.97,0.97,0.84,0.999];
% par = 0.95*ones(1,21);
for i = 1:21
str = num2str(i);
if (i<10)
    str = strcat('0',num2str(i));
end
d01 = load(strcat('D:\研究生\Paper 2\te_process\d',str,'_te.dat'));
% for j = 1:size(d01,1)-wid+1  Dte(j,:) = moving_window(d01,j+wid-1,wid); end
Dte = scale(d01,Dx(1,:),Dx(2,:));
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
control_lim = varx(min(find(fun>par(i))));
Con_lim = control_lim*ones(1,length(Disim_te));
Disim_te = [zeros(1,wid-1),Disim_te];
Con_lim = [zeros(1,wid-1),Con_lim];
int_var = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; fdr(i) = roundn(int_var,-2);
int_var = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  far(i) = roundn(int_var,-2);
D_str{i} = [Disim_te',Con_lim'];
end
multi_plot(D_str);
figure(); hold on 
x = [1:21]; y = [fdr;-far]';
bar(x,y);
sub = legend('Fdr','Far');
set(sub,'FontName','Times','FontSize',13.5);
set(gca,'xtick',[1:22]); plot([x,22],-5*ones(1,22),'k--');plot([x,22],-10*ones(1,22),'k--')
set(gca,'xlim',[0,22]);

%fprintf('fdr = %.2f\n',fdr)

% figure(); hold on 
% plot([1:length(Disim_te)],Disim_te,'b');
% plot([1:length(Disim_te)],control_lim,'r');
% fdr = sum(Disim_te(161:end)>control_lim)/length(Disim_te(161:end))*100; fdr = roundn(fdr,-2);
% far = sum(Disim_te(1:160)>control_lim)/length(Disim_te(1:160))*100;  far = roundn(far,-2);
% title(strcat('FAR:',num2str(far),'%','  FDR:',num2str(fdr),'%'),'FontName','Times','FontSize',13.5);
% figure(); hold on
% stem(h0/sum(h0),'r','fill'); stem(h1/sum(h1),'b','fill');legend('norm','fault')
% axis([0 1000 0 40]);
% set(gca,'yTick',0:10:40);
% lg = test_som(sMap,Dte,V,norm_center)
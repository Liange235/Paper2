function [N_Center,Dx,control_lim] = mysom_train(wid,par)
ind = [1:52];
dnorm = load('D:\ÑĞ¾¿Éú\Paper 2\te_process\d00.dat'); dnorm = dnorm(1:480,:);
for i = 1:size(dnorm,1)-wid+1  D(i,:) = moving_window(dnorm,i+wid-1,wid); end
[D,Dx] = mncn(D);
sDiris = som_data_struct(D,'name','Te (train)');
sMap = som_make(sDiris,'small','tracking',0); 
h0 = som_hits(sMap,sDiris.data);
[~,ind] = sort(h0,'descend');
N_Center = sMap.codebook(ind,:);
sub = sDiris.data-ones(size(D,1),1)*N_Center(1,:);
Disim = diag(sqrt(sub*sub'));
[fun,varx] = ksdensity(Disim,'bandwidth',0.1,'function','cdf');
control_lim = varx(min(find(fun>par)));

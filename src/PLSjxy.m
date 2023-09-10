function [Yq] = PLSjxy(Xtrain,Ytrain,Xtest,Comp)
%% 实现PLS预测
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 准备
[Ntrain,cx]=size(Xtrain);
[~,cy]=size(Ytrain);
[Ntest,~]=size(Xtest);
k=0;                                     %成分数
tol=0.0000000001;

%% 数据标准化
[X,Mx,Sx] = zscore(Xtrain);[Y,My,Sy] = zscore(Ytrain);
Xte = (Xtest - repmat(Mx,Ntest,1))./repmat(Sx,Ntest,1);

%% 参数训练--进行两层迭代
while norm(Y)>tol && k<cx   %外层迭代：确定所有成分，保证个数小于变量数同时能充分解释y
   % 选择随机y的一列作为u
   [dummy,xm] =  max(sum(X.*X));
   [dummy,ym] =  max(sum(Y.*Y));
   u = Y(:,ym);
   t1 = X(:,xm);     %t1用于存放最新的主成分
   t = zeros(Ntrain,1);   %t用于存放旧的主成分
   % 进行内层迭代确定成分
   while norm(t1-t) > tol %内层迭代：确定一个成分，保证收敛完全
       t = t1;
       w = X'*u/((u'*u)); w = w/norm(w);
       t1 = X*w;
       q = Y'*t1/((t1'*t1));
       u = Y*q;
   end
   p = X'*t/(t'*t);

   % 残差阵 
   X = X-t*p';
   Y = Y-t*q';
   % 存放数据
   k = k+1;
   T(:,k) = t1;
   P(:,k) = p;
   W(:,k) = w;
   Q(:,k) = q;
end

%% 数据预测(迭代)
% for K0=1:Comp
%     Yq(:,K0)=zeros(Ntest,1);
%     K1=K0;
%     Xq=Xtest;
%     for i=1:K0
%     t=Xq*W(:,i);
%     Xq=Xq-t*P(:,i)';
%     Yq(:,K1)=Yq(:,K1)+t*B(:,i)';%倒推出Yq
%     end
% end
% Yq = Yq(:,Comp);

%% 数据预测(公式求解)
% Beta = pinv(P(:,1:Comp)*P(:,1:Comp)')*P(:,1:Comp)*diag(B(:,1:Comp))*Q(:,1:Comp)'; %不好
Beta =  W(:,1:Comp)*inv(P(:,1:Comp)'*W(:,1:Comp))*Q(:,1:Comp)';

Yq = Xte*Beta.*repmat(Sy,Ntest,1)+repmat(My,Ntest,1);
end
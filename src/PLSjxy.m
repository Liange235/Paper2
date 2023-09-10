function [Yq] = PLSjxy(Xtrain,Ytrain,Xtest,Comp)
%% ʵ��PLSԤ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ׼��
[Ntrain,cx]=size(Xtrain);
[~,cy]=size(Ytrain);
[Ntest,~]=size(Xtest);
k=0;                                     %�ɷ���
tol=0.0000000001;

%% ���ݱ�׼��
[X,Mx,Sx] = zscore(Xtrain);[Y,My,Sy] = zscore(Ytrain);
Xte = (Xtest - repmat(Mx,Ntest,1))./repmat(Sx,Ntest,1);

%% ����ѵ��--�����������
while norm(Y)>tol && k<cx   %��������ȷ�����гɷ֣���֤����С�ڱ�����ͬʱ�ܳ�ֽ���y
   % ѡ�����y��һ����Ϊu
   [dummy,xm] =  max(sum(X.*X));
   [dummy,ym] =  max(sum(Y.*Y));
   u = Y(:,ym);
   t1 = X(:,xm);     %t1���ڴ�����µ����ɷ�
   t = zeros(Ntrain,1);   %t���ڴ�žɵ����ɷ�
   % �����ڲ����ȷ���ɷ�
   while norm(t1-t) > tol %�ڲ������ȷ��һ���ɷ֣���֤������ȫ
       t = t1;
       w = X'*u/((u'*u)); w = w/norm(w);
       t1 = X*w;
       q = Y'*t1/((t1'*t1));
       u = Y*q;
   end
   p = X'*t/(t'*t);

   % �в��� 
   X = X-t*p';
   Y = Y-t*q';
   % �������
   k = k+1;
   T(:,k) = t1;
   P(:,k) = p;
   W(:,k) = w;
   Q(:,k) = q;
end

%% ����Ԥ��(����)
% for K0=1:Comp
%     Yq(:,K0)=zeros(Ntest,1);
%     K1=K0;
%     Xq=Xtest;
%     for i=1:K0
%     t=Xq*W(:,i);
%     Xq=Xq-t*P(:,i)';
%     Yq(:,K1)=Yq(:,K1)+t*B(:,i)';%���Ƴ�Yq
%     end
% end
% Yq = Yq(:,Comp);

%% ����Ԥ��(��ʽ���)
% Beta = pinv(P(:,1:Comp)*P(:,1:Comp)')*P(:,1:Comp)*diag(B(:,1:Comp))*Q(:,1:Comp)'; %����
Beta =  W(:,1:Comp)*inv(P(:,1:Comp)'*W(:,1:Comp))*Q(:,1:Comp)';

Yq = Xte*Beta.*repmat(Sy,Ntest,1)+repmat(My,Ntest,1);
end
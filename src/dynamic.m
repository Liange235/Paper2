function y = dynamic(x,t)
% ����x��ʱ�Ӵ���t���Ծ���x������ʱ�Ӷ�̬��
% ���y
% ����x belongs to mxn
% 1.�����趨ͳһ��ʱ������t��x���ж�̬������y��y belongs to (m-t)x[n(t+1)]
%  
[m,n] = size(x);
y  = [];
temp = 1:m-t;
lab = 1:m-t;
for i = 1:t
     temp = temp+ones(1,m-t);
     lab = [lab;temp];
end
lab  = lab';
for i = 1:n
    temp = x(:,i);
    y = [y temp(lab)];
end
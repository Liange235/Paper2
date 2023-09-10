%%%%%%%%%******moving window approach******%%%%%%%%%
function y = moving_window(x,n1,n2)
if nargin<3
    start = size(x,1);
end
if n1<n2
    error('not enough data in the window');
end
start = n1; data = x; width = n2;
y = [mean(data(start-width+1:start,:))];
if n2==1
    y = [data(start-width+1:start,:)];
end


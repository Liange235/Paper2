%%**plot a picture involving different clusters according to the labels
function [] = myscatter(data,label)
cluster_no = max(label);
color = [255,0,0;0,255,0;0,0,255;255,222,0;248,147,29;119,52,96;145,145,141;160,191,124;36,169,225;89,69,61;230,16,155]*1/255;
type = 'opx*s^+dh.>';
% type = 'odh';
hold on
for i = 1:cluster_no
    m = type(i);
    ind = find(label == i);
    c = repmat(color(i,:),length(ind),1);
    scatter(data(ind,1),data(ind,2),140,c,m,'LineWidth',0.5);
end
%%%**************test color code
% label = [1:11]';
% data = randn(11,2);
% data = bsxfun(@times,data,label);
% myscatter(data,label);
% legend('1','2','3','4','5','6','7','8','9','10','11');
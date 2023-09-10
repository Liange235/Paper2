function [] = multi_plot(input)
p = length(input); 
sub = factor(p);
m = sub(1);n = sub(2);
if(length(factor(p))>2)
    for i = 1:length(factor(p))
        f1 = prod(sub(1:i));
        f2 = p/f1;
        if f1>f2
            break
        end
    end
m = prod(sub(1:i-1)); n = p/m;
m = f1*(abs(f1-f2)<abs(m-n))+m*(abs(f1-f2)>=abs(m-n)); n = p/m; 
end
figure()
for i = 1:p
   Disim_te = input{i};
   control_lim = max(Disim_te(:,2));
   Disim_te = Disim_te(:,1);
   s = subplot(m,n,i); hold on
   plot([1:length(Disim_te)],Disim_te,'b');
   plot([1:length(Disim_te)],control_lim,'r');
   title(s,strcat('Fault',num2str(i)),'FontName','Times','FontSize',13.5);
   plot([160,160],[0,max(Disim_te)+10],'k--')
end

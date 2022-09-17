clc
clear
fileID = fopen('count1.bin');
count1 = fread(fileID,[1,101],'double')

fileID = fopen('Prob.out');
Pvec = fread(fileID,[1,101],'double')

% plot(count1); 
figure (1)
plot(Pvec,'k','linewidth',2)
xlabel('k'); ylabel('p(k/10)');
title('Probability')
xlim([0,length(Pvec)])
set(gca,'fontsize',20)
clc
clear
fileID = fopen('Xvec.out');
Xvec2 = fread(fileID,[1,1001],'double')

fileID = fopen('Yvec.bin');
Yvec2 = fread(fileID,[1,1001],'double')

figure (1)
hold on
plot(Xvec2,Yvec2,'r','linewidth',2); hold on
xlabel('x');
ylabel('y');
title('Phase plane plots')
set(gca,'fontsize',20)
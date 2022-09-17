clc
clear
fileID = fopen('finalu.bin');
sol = fread(fileID,[1,100],'double');
figure (1)
plot(sol)

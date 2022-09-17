clc
clear
N = 128;
M = 2*N-1;
count=5000;
fileID = fopen('sources.out');
sol = fread(fileID,[M-2,N-2],'double');

xvec = linspace(-2,2,M-2);
yvec = linspace(-1,1,N-2);
[X,Y] = meshgrid(xvec,yvec);
size(X)
size(sol)
figure (1)
mesh(X,Y,sol');
xlabel('x')
ylabel('y')
title('Numeircal solution with N = 128 and \lambda = 100')
set(gca,'fontsize',20)

% figure (2)
% fileID = fopen('resmat.bin');
% res = fread(fileID,[M-2,N-2],'double');
% mesh(res)
% norm(sol)
% norm(res)

% figure (2)
% fileID = fopen('residual.bin');
% resnorm = fread(fileID,[count,1],'double');
% plot(resnorm)


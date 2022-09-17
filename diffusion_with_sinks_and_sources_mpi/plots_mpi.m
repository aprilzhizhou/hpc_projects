clc
clear
N = 128;
M = 2*N-1;
Ns = N/2+1;
count=10000;
fileID = fopen('finalu.bin');
sol = fread(fileID,[1,N*M],'double');
sol = reshape(sol,[M,N]);
% xvec = linspace(-2,2,M-2);
% yvec = linspace(-1,1,N-2);
% [X,Y] = meshgrid(xveMc,yvec);
% size(X)
% size(sol)
figure (1)
surf(sol,'EdgeColor','none');
xlabel('x')
ylabel('y')
% title('Numeircal solution with N = 128 and \lambda = 100')
% set(gca,'fontsize',20)

% figure (2)
% fileID = fopen('resmat.bin');
% res = fread(fileID,[M-2,N-2],'double');
% mesh(res)
% norm(sol)
% norm(res)

figure (2)
fileID = fopen('residual.out');
resmax = fread(fileID,[count,1],'double');
plot(resmax);
resmax(end)

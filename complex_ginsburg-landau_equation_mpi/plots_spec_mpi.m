clc
clear
N = 128;
 
fileID2 = fopen('CGL.out');
solt = fread(fileID2,[1,N*N*2*10],'double');
soltreal= solt(1:2:end);
soltimag = solt(2:2:end);
soltreal = reshape(soltreal,[N,N*10]);
soltimag = reshape(soltimag,[N,N*10]);
contourf(soltreal,'EdgeColor','none');

% 
% 
% fileID1 = fopen('finalsol.bin');
% sol = fread(fileID1,[1,N*N*2],'double');
% solreal= sol(1:2:end);
% solimag = sol(2:2:end);
% solreal = reshape(solreal,[N,N]);
% solimag = reshape(solimag,[N,N]);
% mag = sqrt(solreal.^2 + solimag.^2);
% 
% 
% % xvec = linspace(0,128*pi,N);
% % yvec = linspace(0,128*pi,N);
% % [X,Y] = meshgrid(xvec,yvec);
% 
% figure (1)
% subplot(1,2,1)
% contourf(solreal,'EdgeColor','none');
% title("Re(A)")
% set(gca,'fontsize',18)
% axis equal
% subplot(1,2,2)
% contourf(solimag,'EdgeColor','none');
% title("Im(A)")
% set(gca,'fontsize',18)
% axis equal
% figure (2)
% subplot(1,2,1)
% contourf(mag,'EdgeColor','none');
% title("|A|")
% set(gca,'fontsize',18)
% axis equal
% subplot(1,2,2)
% phase = atan2(solimag,solreal);
% contourf(phase,'EdgeColor','none');
% title("arctan(Im(A)/Re(A))")
% set(gca,'fontsize',18)
% axis equal
% 
% 
% 

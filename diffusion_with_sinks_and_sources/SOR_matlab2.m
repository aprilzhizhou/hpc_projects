%% initialization
w = 1.8;
N =128;
M = 2*N-1;
Uinit = zeros(M,N);
tol = 1e-4;
lam = 10;
kMax = 1000;
f = @(x,y) 10*lam/sqrt(pi)*exp(-lam^2*((x-1).^2+y.^2)) ...
     - 10*lam/sqrt(pi)*exp(-lam^2*((x+1).^2+y.^2));

% construct the solution matrix
U = zeros(size(Uinit));

h = 2/(N-1);
U = Uinit; Res = [];
tic
for n = 1:kMax
    for jj = 1
        for kk = 1
            U(jj,kk)  = 0; 
            
        end
        for kk = 2:N-1
            U(jj,kk)  = U(jj+1,kk);
        end
        for kk = N
            U(jj,kk)  = 0;
        end
    end
    
    for jj = 2:M-1
        for kk = 1
            U(jj,kk)  = 0;
        end
        for kk = 2:N-1
            U(jj,kk)  = (1-w)*U(jj,kk)               ...
                +w/4*(           U(jj,kk+1)                 ...
                + U(jj-1,kk)                + U(jj+1,kk) ...
                +U(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
            res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                  +U(jj,kk-1)+ U(jj-1,kk)+U(jj+1,kk)+U(jj,kk+1));
        end
        for kk = N
            U(jj,kk)  = 0;
        end
    end
    
    
    for jj = M
        for kk = 1
            U(jj,kk)  = 0;
        end
        for kk = 2:N-1
            U(jj,kk)  = U(jj-1,kk);
        end
        for kk = N
            U(jj,kk)  = 0;
        end
    end
   
    Res = [Res; norm(res,inf)];
        if Res(end) < tol
        break
    end
  
end
toc 
figure (1)
mesh(U(1:end,1:end));
% Res
norm(U)
norm(res)
figure (2)
plot(Res)
%     for ii = 1:M-2
% %         for jj = 2:N+1
% %             U(ii,jj) = w/(-4+h^2*pi^2*K)*(h^2*F(ii,jj) ...
% %                 -(U(ii,jj-1)+ U(ii-1,jj)+U(ii+1,jj)+U(ii,jj+1)))...
% %         + (1-w)*U(ii,jj);
%         end
%     end
%     for ii = M-1
%     end
%     for ii = 2:N+1
%         for jj = 2:N+1
%             res(ii,jj) = F(ii,jj) -1/h^2*( (-4 + h^2*pi^2*K)*U(ii,jj) ...
%                 +U(ii,jj-1)+ U(ii-1,jj)+U(ii+1,jj)+U(ii,jj+1));
%         end
%     end
%     Res = [Res; norm(res,inf)];
%     if Res(end) < tol
%         break
%     end
%     


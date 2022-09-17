%% initialization
w = 1.8;
N =128;
M = 2*N-1;
Uinit = zeros(M-1,N-1);
tol = 1e-9;
lam = 10;
kMax = 1000;
 f = @(x,y) 10*lam/sqrt(pi)*exp(-lam^2*((x-1).^2+y.^2)) ...
     - 10*lam/sqrt(pi)*exp(-lam^2*((x+1).^2+y.^2));
% f = @(x,y) exp(-1/2*x-1/2*y);
% construct the solution matrix
Unew = zeros(size(Uinit));

h = 2/(N-1);
Uold = Uinit; Res = [];
res = zeros(size(Uinit));
tic
for n = 1:kMax
    for jj = 1
        for kk = 1
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                + Unew(jj,kk)                + Uold(jj+1,kk) ...
                +0)                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = 2:N-3
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                + Unew(jj,kk)                + Uold(jj+1,kk) ...
                +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = N-2
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           0                 ...
                + Unew(jj,kk)                + Uold(jj+1,kk) ...
                +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
    end
    
    for jj = 2:M-2
        for kk = 1
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                + Unew(jj-1,kk)                + Uold(jj+1,kk) ...
                                 +0)                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = 2:N-3
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                + Unew(jj-1,kk)                + Uold(jj+1,kk) ...
                +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = N-2
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           0                 ...
                + Unew(jj-1,kk)                + Uold(jj+1,kk) ...
                +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
    end
    
    
    for jj = M-2
        for kk = 1
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                + Unew(jj-1,kk)                + Uold(jj,kk) ...
                +0)                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = 2:N-3
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           Uold(jj,kk+1)                 ...
                +Unew(jj-1,kk)                + Uold(jj,kk) ...
                +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
        for kk = N-2
            Unew(jj,kk)  = (1-w)*Uold(jj,kk)               ...
                +w/4*(           0                 ...
                + Unew(jj-1,kk)                + Uold(jj,kk) ...
                              +Unew(jj,kk-1))                ...
                - w/4*h^2*f(-2+jj*h,-1+kk*h);
        end
    end
    % compute residuals 
    for jj = 1
        for kk = 1
     res(jj,kk) = -h^2/4 * f(-2+jj*h,-1+kk*h) +1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk+1)+ Unew(jj,kk)+Unew(jj+1,kk));
        end
        for kk = 2:N-3
            res(jj,kk) = -h^2/4 * f(-2+jj*h,-1+kk*h) +1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj,kk)+Unew(jj+1,kk)+Unew(jj,kk+1));
        end
        for kk = N-2
        res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj,kk)+Unew(jj+1,kk)+0);
        end
    end
    
    for jj = 2:M-2
        for kk = 1
            res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +0+ Unew(jj-1,kk)+Unew(jj+1,kk)+Unew(jj,kk+1));
        end
        for kk = 2:N-3
           res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj-1,kk)+Unew(jj+1,kk)+Unew(jj,kk+1));
        end
        for kk = N-2
        res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj-1,kk)+Unew(jj+1,kk)+0);
        end
    end
    
    
    for jj = M-2
        for kk = 1
           res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +0+ Unew(jj-1,kk)+Unew(jj,kk)+Unew(jj,kk+1));
        end
        for kk = 2:N-3
    res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj-1,kk)+Unew(jj,kk)+Unew(jj,kk+1));
        end
        for kk = N-2
    res(jj,kk) = -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*Unew(jj,kk) ...
                  +Unew(jj,kk-1)+ Unew(jj-1,kk)+Unew(jj,kk));
        end
    end
    
    
    
%     for jj = 2:M-3
%         for kk = 2:N-2
%            res(jj,kk) = f(-2+jj*h,-1+kk*h) -1/h^2*(Unew(jj,kk) ...
%                  +Unew(jj,kk-1)+ Unew(jj-1,kk)+Unew(jj+1,kk)+Unew(jj,kk+1));
%         end
%     end
    Res = [Res; norm(res,inf)];
        if Res(end) < tol
        break
    end
   
    Uold = Unew;
end
toc 
figure (1)
mesh(Unew(1:end-2,1:end-2));
% Res
norm(Unew)
norm(res)
figure (2)
plot(Res)
%     for ii = 1:M-2
% %         for jj = 2:N+1
% %             Unew(ii,jj) = w/(-4+h^2*pi^2*K)*(h^2*F(ii,jj) ...
% %                 -(Unew(ii,jj-1)+ Unew(ii-1,jj)+Uold(ii+1,jj)+Uold(ii,jj+1)))...
% %         + (1-w)*Uold(ii,jj);
%         end
%     end
%     for ii = M-1
%     end
%     for ii = 2:N+1
%         for jj = 2:N+1
%             res(ii,jj) = F(ii,jj) -1/h^2*( (-4 + h^2*pi^2*K)*Unew(ii,jj) ...
%                 +Unew(ii,jj-1)+ Unew(ii-1,jj)+Unew(ii+1,jj)+Unew(ii,jj+1));
%         end
%     end
%     Res = [Res; norm(res,inf)];
%     if Res(end) < tol
%         break
%     end
%     


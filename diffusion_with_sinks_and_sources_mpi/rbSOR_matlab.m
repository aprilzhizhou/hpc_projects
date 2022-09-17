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
res = zeros(size(U));
tic


for n = 1:kMax
    %% update red nodes
    for jj =1 % BC: du/dx(-2,y) =0;
        for kk = 2:N-1
            if (mod(jj+kk,2)==1) % update red values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj,kk)...
                    + U(jj+1,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj,kk)+U(jj+1,kk)+U(jj,kk+1));
            end
        end
    end
    for jj=2:M-1
        for kk = 2:N-1
            if (mod(jj+kk,2)==1) % update red values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj-1,kk)...
                    + U(jj+1,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj-1,kk)+U(jj+1,kk)+U(jj,kk+1));
            end
        end
    end
    for jj = M
        for kk = 2:N-1
            if (mod(jj+kk,2)==1) % update red values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj-1,kk)...
                    + U(jj,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj-1,kk)+U(jj,kk)+U(jj,kk+1));
            end
        end
    end
    
    %% update black nodes 
        for jj =1 % BC: du/dx(2,y) =0;
        for kk = 2:N-1
            if (mod(jj+kk,2)==0) % update black values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj,kk)...
                    + U(jj+1,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj,kk)+U(jj+1,kk)+U(jj,kk+1));
            end
        end
    end
    for jj=2:M-1
        for kk = 2:N-1
            if (mod(jj+kk,2)==0) % update black values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj-1,kk)...
                    + U(jj+1,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj-1,kk)+U(jj+1,kk)+U(jj,kk+1));
            end
        end
    end
        for jj = M
        for kk = 2:N-1
            if (mod(jj+kk,2)==0) % update black values
                U(jj,kk)  = (1-w)*U(jj,kk)+w/4*(U(jj,kk+1)+ U(jj-1,kk)...
                    + U(jj,kk)+U(jj,kk-1))- w/4*h^2*f(-2+jj*h,-1+kk*h);
                res(jj,kk) =  -h^2/4*f(-2+jj*h,-1+kk*h) + 1/4*(-4*U(jj,kk) ...
                    +U(jj,kk-1)+ U(jj-1,kk)+U(jj,kk)+U(jj,kk+1));
            end
        end
    end
    Res = [Res max(max(res))];
end

toc
figure (1)
mesh(U(1:end,1:end));

% norm(U)
% norm(res)
figure (2)
plot(Res)
Res(end)


N = 32;
A = ones(N,N);
Alap = zeros(N,N);
tic
for i = 1:N
    for j = 1:N
        ii = i-1;
        jj = j-1;
        if (ii <N/2 && jj <N/2)
            Alap(i,j) =-(ii^2 + jj^2)*A(i,j) ;
        end
        if (ii < N/2 && jj >= N/2)
            Alap(i,j) =-(ii^2 + (jj-N)^2)*A(i,j) ;
        end
        if (ii >= N/2 && jj < N/2)
            Alap(i,j) =-((ii-N)^2 + jj^2)*A(i,j) ;
        end
        if (ii >= N/2 && jj >= N/2)
            Alap(i,j) =-((ii-N)^2 + (jj-N)^2)*A(i,j) ;
        end
    end
end
toc

tic
Alap1 = zeros(N,N);
for i = 1:N
    for j = 1:N
        ii = i-1;
        jj = j-1;
        if jj >= N/2 
            jj1 = jj-N;
        end
        if ii >= N/2
            ii1 = ii-N;
        end
        Alap1(i,j) =-(ii1^2 + jj1^2)*A(i,j) ;
    end
end
toc 

norm(Alap-Alap1)
alpha =1; sigma = 1/2; N = 10000; M = 1000;
T = 10; dt = T/M; K = M/10 + 1;
count = zeros(M/10+1,1);

for n =0:N-1
   u1 = rand(); u2 = rand();
   X = sqrt(-2*log(u1))*cos(2*pi*u2);
   Y =0;
   if norm([X+alpha,Y]) <= alpha/2 || norm([X-alpha,Y]) <= alpha/2
       count(1) = count(1)+1;
   end
   
   for m = 0:M-1
       t = (m+1)*dt;
       v1 = rand(); v2 = rand();
       dW = dt*sqrt(-log(v1))*cos(2*pi*v2);
       X = X + dt*Y;
       Y = Y + ((alpha^2-X^2)*X-Y)*dt + sigma*X*dW;
       if mod(m+1,10) ==0
           if norm([X+alpha,Y]) <= alpha/2 || norm([X-alpha,Y]) <= alpha/2
               l = (m+1)/10 + 1;
               count(l) = count(l)+1;
           end
       end
   end
end
count 
plot(count)
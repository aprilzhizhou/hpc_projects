N = [25000, 50000, 100000, 200000]';
T = [1.610845
3.426106
6.073540
11.258840];


figure (1)
loglog(N,T,'b-o');
hold on
loglog(N,N,'g-+')
loglog(N,N.^2,'r-x')
legend({'Computatin time','N','N^2'},'location','best')
set(gca,'fontsize',18)

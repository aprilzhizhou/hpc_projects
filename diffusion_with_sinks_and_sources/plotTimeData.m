t64 = 0.000417;
t128=0.001786;
t256 = 0.009424;
figure (1)
T=[t64 t128 t256];
N=[64,128,256];
loglog(N,T,'b-o');
hold on
loglog(N,N,'g-+')
loglog(N,N.^2,'r-x')
legend({'Timing data','N','N^2'},'location','best')
set(gca,'fontsize',18)

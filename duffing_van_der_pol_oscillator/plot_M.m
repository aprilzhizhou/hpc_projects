M = [250, 500, 1000, 2000]';
T = [1.375889
2.722083
5.403277
10.671555];


figure (1)
loglog(M,T,'b-o');
hold on
loglog(M,M,'g-+')
loglog(M,M.^2,'r-x')
legend({'Computatin time','M','M^2'},'location','best')
set(gca,'fontsize',18)

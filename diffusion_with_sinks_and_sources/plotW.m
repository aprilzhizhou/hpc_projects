wvec = [1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98]';
itervec = [905 
821 
731 
644 
560 
477 
392 
294 
371 
706  ];
figure (1)
plot(wvec, itervec,'k','linewidth',2);
title('Iter. nubmer vs. \omega with N = 128');
set(gca,'fontsize',18)
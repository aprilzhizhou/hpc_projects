This projects uses the red-black SOR method (with relaxation parameter $\omega =1.8$) to solve the elliptic equation 
$\nabla^2 u = \frac{10\lambda}{\sqrt{\pi}} (\exp{(-\lambda^2 ((x-1)^2+y^2) )} -\exp{(-\lambda^2 ((x+1)^2+y^2) )} ),$
with boundary conditions 
$$\frac{\partial u}{\partial x}(\pm2,y) = 0,\quad u(x,\pm 1) = 0,$$
for $-2 \leq x \leq 2, -1 \leq y \leq 1,$ and $\lambda = 100.$ The numerical method uses a $(2N-1)\times N$ grid with red-black ordering and `tol = 1e-9` as the stopping criterion. 
The numerical implementation in this program is parallel using the MPI module.

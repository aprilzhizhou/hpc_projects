This project uses the pseudospetral method with fourth-order Runge-Kutta mehtod to solve the complex Ginsburg-Landau equation, which is given by 
$$\frac{\partial A}{\partial t} = A + (1+ic_1)\nabla^2 A - (1-c_3)|A|^2 A,$$
where $A(x,t)$ is a complex valued field. 
The program solves the complex Ginsburg-Landau equation in 2D on rectangular domain with length $L = 128\pi$ on each side. Converting to a length of $2\pi$ on each side,
the equation can be rescaled to be 
$$\frac{\partial A}{\partial t} = A + (\frac{2\pi}{L})^2 (1+ic_1)\nabla^2 A - (1-c_3)|A|^2 A.$$

We use the psedospetral method with RK4 with random initial data in the range $[-1.5,1.5] + i[-1.5,1.5]$ and run until terminal time $T= 10^4$. For parameters, $c_1  = 1.5, 0< c_3<0.75.$

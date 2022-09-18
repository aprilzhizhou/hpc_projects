This project uses the pseudospetral method with fourth-order Runge-Kutta mehtod to solve the complex Ginsburg-Landau equation, which is given by 
$$\frac{\partial A}{\partial t} = A + (1+ic_1)\nabla^2 A - (1-c_3)|A|^2 A,$$
where $A(x,t)$ is a complex valued field. 
The program solves the complex Ginsburg-Landau equation in 2D on rectangular domain with length $L = 128\pi$ on each side. Converting to a length of $2\pi$ on each side,
the equation can be rescaled to be 
$$\frac{\partial A}{\partial t} = A + (\frac{2\pi}{L})^2 (1+ic_1)\nabla^2 A - (1-c_3)|A|^2 A.$$

We use the psedospetral method with RK4 with random initial data in the range $[-1.5,1.5] + i[-1.5,1.5]$ and run until terminal time $T= 10^4$. For parameters, $c_1  = 1.5, 0< c_3<0.75,$ spiral wave pattern should emerge.  

The program takes four arguments with an optional fifth argument for the random seed so that the program `cgl` can be executed like `$ cgl 128 1.5 0.25 100000`,
where $N=128$ is the size of the domain in $x$ and $y$ directions, $c_1 = 1.5$ and $c_3 = 0.25$ arethe real-valued parameters, and $M = 100,000$ is the number of time steps. The optional fifth argument is an integer seed value for generating random initial values. The program prints the grid values of $A$ in a file called `CGL.out` that contains the data at the time points $t = 1000k$ for $k = 0, 1, \cdots, 10.$

The numerical method uses an $N\times N$ grild. The implemential is parallel with the MPI module. 

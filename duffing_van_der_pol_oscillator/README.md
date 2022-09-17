
This project uses the Euler-Maruyama method to solve the Duffing-Van der Pol equation, which can be written as a first-order system of stochastic differential equations given by
$$dX_t = Y_t dt, \quad dY_t ((\alpha^2 - X_t^2)X_t - Y_t) + \sigma X_t dW_t. $$
The program takes four arguments plus an optional argument for the seed so that the program named vdp can be run with the commend`$ vdp 1 0.5 100 100000`, 
where $\alpha$ = 1 is the half distance between the two stable attractors, $\sigma = 0.5$ is the strength of the noise, $M=1000$ is the number of time steps to use, and $N = 100,000$ is the number of trials. The terminal time for each trial is $T = 10$ so that the time step size is $\Delta t = T/M.$ Define the function $p(t)$ to be the probability of being close to an equilibrium point by 
$$p(t) = P (  ||(X_t,Y_t) -(-\alpha,0)||  \leq \frac{\alpha}{2})  +  P ( ||(X_t,Y_t) -(\alpha,0)|| \leq \frac{\alpha}{2} ).$$
For each time value $t = k/10$, the program counts how many trials resulted in the point $(X_t,Y_t)$ satisfying one of the two above inequalities and creates a double precision array `P[101]` where `P[k]` $\approx p(k/10)$ for $0 \leq k \leq 100.$ The array `P` is saved to a file named "`Prob.out`" in binary format.  

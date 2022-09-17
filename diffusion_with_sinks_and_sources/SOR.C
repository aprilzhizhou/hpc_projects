#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double fun(double lam, double x, double y);
double fun(double lam, double x, double y)
{
    return  10*lam/sqrt(M_PI)*exp(-lam*lam*((x-1)*(x-1)+y*y))-10*lam/sqrt(M_PI)*exp(-lam*lam*((x+1)*(x+1)+y*y));
    // return exp(-0.5*x -0.5* y);
}



int main(int argc, const char * argv[])
{
    double w = 1.8;
    int N =201;
    int M = 2*N-1;
    double Uold[M-2][N-2];
    double Unew[M-2][N-2];
    double tol = 1e-6;
    double lam = 3;
    double h = 2./200.;
    printf("h = %10.10f",h);
    int K = 100;
    
    double fvec[N-2][M-2];
    for (int m=0; m<=M-3; ++m)
    {
        for (int n=0; n<=N-3; ++n)
        {
            fvec[n][m] = fun(lam,-2+(m+1)*h,-1+(n+1)*h);
        }
    }
    FILE *funct = fopen("funct.bin","w");
    fwrite(fvec,sizeof(double),(M-2)*(N-2),funct);
    fclose(funct);
    
    //printf("h = %10.10f \n ",h);
    // initialize Uold
    for (int m=0; m<M-2; ++m)
    {
        for (int n=0; n<N-2; ++n)
        {
            Uold[m][n] = 0;
        }
    }
    
    
    for (int count=0; count<K; ++count)
    {
        //j = 0;
        //k = 0;
        Unew[0][0]  = (1-w)*Uold[0][0]
        +w/4*(                Uold[0][1]
              + Unew[0][0]                + Uold[1][0]
              +    0       )
        - w/4*h*h*fun(lam,-2+(1)*h,-1+(1)*h);
        for (int k=1; k<=N-4; ++k)
        {
            Unew[0][k]  = (1-w)*Uold[0][k]
            +w/4*(                Uold[0][k+1]
                  + Unew[0][k]                + Uold[1][k]
                  +Unew[0][k-1])
            - w/4*h*h*fun(lam,-2+(1)*h,-1+(k+1)*h);
        }
        //k=N-3;
        Unew[0][N-3]  = (1-w)*Uold[0][N-3]
        +w/4*(                    0
              + Unew[0][N-3]                + Uold[1][N-3]
              +Unew[0][N-4])
        - w/4*h*h*fun(lam,-2+(1)*h,-1+(N-2)*h);
        for (int j=1; j<=M-4; ++j)
        {
            //k=0;
            Unew[j][0]  = (1-w)*Uold[j][0]
            +w/4*(                Uold[j][1]
                  + Unew[j-1][0]                + Uold[j+1][0]
                  +    0       )
            - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(1)*h);
            for (int k=1; k<=N-4; ++k)
            {
                Unew[j][k]  = (1-w)*Uold[j][k]
                +w/4*(                Uold[j][k+1]
                      + Unew[j-1][k]                + Uold[j+1][k]
                      +Unew[j][k-1])
                - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(k+1)*h);
            }
            // k=N-3;
            Unew[j][N-3]  = (1-w)*Uold[j][N-3]
            +w/4*(                    0
                  + Unew[j-1][N-3]                + Uold[j+1][N-3]
                  +Unew[j][N-4])
            - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(N-2)*h);
        }
        // j = M-3;
        // k = 0;
        Unew[M-3][0]  = (1-w)*Uold[M-3][0]
        +w/4*(                Uold[M-3][1]
              + Unew[M-4][0]                + Uold[M-3][0]
              +    0       )
        - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(1)*h);
        for (int k=1; k<=N-4; ++k)
        {
            Unew[M-3][k]  = (1-w)*Uold[M-3][k]
            +w/4*(                Uold[M-3][k+1]
                  + Unew[M-4][k]                + Uold[M-3][k]
                  +Unew[M-3][k-1])
            - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(k+1)*h);
        }
        // k=N-3;
        Unew[M-3][N-3]  = (1-w)*Uold[M-3][N-3]
        +w/4*(                    0
              + Unew[M-4][N-3]                + Uold[M-3][N-3]
              +Unew[M-3][N-4])
        - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(N-2)*h);
        
        for (int p=0; p<=M-3; ++p)
        {
            for (int q=0; q<=N-3; ++q)
            {
                Uold[p][q] = Unew[p][q];
            }
        }
    }
    printf("Count = %d",count)
    // write file
//    FILE *sol = fopen("sol.bin","w");
//    fwrite(Unew,sizeof(double),(M-2)*(N-2),sol);
//    fclose(sol);
    return 0;
}





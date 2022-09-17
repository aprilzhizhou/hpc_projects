
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

double fun(double lam, double x, double y);
double findNorm(int M, int N, double res[M][N] );

double fun(double lam, double x, double y)
{
    return  10*lam/sqrt(M_PI)*exp(-lam*lam*((x-1)*(x-1)+y*y))-10*lam/sqrt(M_PI)*exp(-lam*lam*((x+1)*(x+1)+y*y));
}


double findNorm(int M, int N, double res[N][M])
{
    double sum = 0.;
    for (int m=0; m<M; ++m)
    {
        for (int n=0; n<N; ++n)
        {
            sum = sum + res[n][m]*res[n][m];
        }
    }
    return sqrt(sum);
}


int main(int argc, const char * argv[])
{
    double w = 1.8;
    int N =128;
    int M = 2*N-1;
    
    double* U = (double*)malloc(M*N * sizeof(double));
    double res[M][N];
    double tol = 1e-9;
    double lam = 10;
    double h = 2./(N-1);
    printf("h = %10.10f\n",h);
    int K = 2000;
    int count = 0;
    double* resmaxvec = (double*)malloc(K*sizeof(double));
    double resmax = 1000.;
    double val;
    
    // calculate time - initialize time variables
    struct timeval t_start;
    struct timeval t_end;
    
    // initialize U
    for (int j=0; j<=M-1; ++j)
    {
        for (int k=0; k<=N-1; ++k)
        {
            U[N*j+k] = 0;
        }
    }
    
    // initialize res
    for (int j=0; j<=M-1; ++j)
    {
        for (int k=0; k<=N-1; ++k)
        {
            res[j][k] = 0;
        }
    }
    
    //    double* F = (double*)malloc(M*N * sizeof(double));
    //    for (int j=0; j<M; ++j)
    //    {
    //        for (int k=0; k<N; ++k)
    //        {
    //            F[N*j+k] = fun(lam,-2+(j+1)*h,-1+(k+1)*h);
    //            printf("value = %f\n",-2+j*h);
    //        }
    //    }
    //    FILE *Fun = fopen("F.out","w");
    //    fwrite(F,sizeof(double),M*N,Fun);
    //    fclose(Fun);
    
    // get starting time
    gettimeofday(&t_start, 0);
    
    while (count < K && resmax > tol)
    {
        resmax = 0;
        // update red nodes
        // j = 0
        for (int k=1; k<=N-2; ++k)
        {
            if ((0+k)%2 ==0 )
            {
                U[N*0+k] = (1-w)*U[N*0+k]+0.25*w*(U[N*0+(k+1)]+ U[N*(0)+k]+ U[N*(0+1)+k]+U[N*0+(k-1)] )- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
            }
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k=1; k<=N-2; ++k)
            {
                if ((j+k)%2 ==0 )
                {
                    val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[N*j+k]+U[N*j+(k-1)]+U[N*(j-1)+k]+ U[N*(j+1)+k] + U[N*j+(k+1)]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    U[N*j+k] = (1-w)*U[N*j+k]+0.25*w*(U[N*j+(k+1)]+ U[N*(j-1)+k]+ U[N*(j+1)+k]+U[N*j+(k-1)] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
                }
            }
        }
        // j=M-1
        for (int k=1; k<=N-2; ++k)
        {
            if ((M-1+k)%2 ==0 )
            {
                U[N*(M-1)+k] = (1-w)*U[N*(M-1)+k]+0.25*w*(U[N*(M-1)+(k+1)]+ U[N*(M-1-1)+k]+ U[N*(M-1)+k]+U[N*(M-1)+(k-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
            }
        }
        // update black nodes
        // j = 0
        for (int k=1; k<=N-2; ++k)
        {
            if ((0+k)%2 == 1 )
            {
                U[N*0+k] = (1-w)*U[N*0+k]+0.25*w*(U[N*0+(k+1)]+ U[N*(0)+k]+ U[N*(0+1)+k]+U[N*0+(k-1)] )- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
            }
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k=1; k<=N-2; ++k)
            {
                if ((j+k)%2 == 1 )
                {
                    U[N*j+k] = (1-w)*U[N*j+k]+0.25*w*(U[N*j+(k+1)]+ U[N*(j-1)+k]+ U[N*(j+1)+k]+U[N*j+(k-1)] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
                    val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[N*j+k]+U[N*j+(k-1)]+U[N*(j-1)+k]+ U[N*(j+1)+k] + U[N*j+(k+1)]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                }
            }
        }
        // j=M-1
        for (int k=1; k<=N-2; ++k)
        {
            if ((M-1+k)%2 ==1 )
            {
                U[N*(M-1)+k] = (1-w)*U[N*(M-1)+k]+0.25*w*(U[N*(M-1)+(k+1)]+ U[N*(M-1-1)+k]+ U[N*(M-1)+k]+U[N*(M-1)+(k-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
            }
        }
        // compute norm of resisual
        //        resnorm =findNorm(&res[0][0],M-2,N-2);
//        resmax = findNorm(M,N,res);
        resmaxvec[count] = resmax;
        count++;
    }
    
    
    //    while ( count < K && resnorm > tol)
    //    {
    //        //j = 0;
    //        //k = 0;
    //        Unew[0][0]  = (1-w)*Uold[0][0]
    //        +w/4*(                Uold[1][0]
    //              + Unew[0][0]                + Uold[0][1]
    //              +    0       )
    //        - w/4*h*h*fun(lam,-2+(1)*h,-1+(1)*h);
    //        for (int k=1; k<=N-4; ++k)
    //        {
    //            Unew[k][0]  = (1-w)*Uold[k][0]
    //            +w/4*(                Uold[k+1][0]
    //                  + Unew[k][0]                + Uold[k][1]
    //                  +Unew[k-1][0])
    //            - w/4*h*h*fun(lam,-2+(1)*h,-1+(k+1)*h);
    //        }
    //        //k=N-3;
    //        Unew[N-3][0]  = (1-w)*Uold[N-3][0]
    //        +w/4*(                    0
    //              + Unew[N-3][0]              + Uold[N-3][1]
    //              +Unew[N-4][0])
    //        - w/4*h*h*fun(lam,-2+(1)*h,-1+(N-2)*h);
    //        for (int j=1; j<=M-4; ++j)
    //        {
    //            //k=0;
    //            Unew[0][j]  = (1-w)*Uold[0][j]
    //            +w/4*(                Uold[1][j]
    //                  + Unew[0][j-1]                + Uold[0][j+1]
    //                  +    0       )
    //            - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(1)*h);
    //            for (int k=1; k<=N-4; ++k)
    //            {
    //                Unew[k][j]  = (1-w)*Uold[k][j]
    //                +w/4*(                Uold[k+1][j]
    //                      + Unew[k][j-1]                + Uold[k][j+1]
    //                      +Unew[k-1][j])
    //                - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(k+1)*h);
    //            }
    //            // k=N-3;
    //            Unew[N-3][j]  = (1-w)*Uold[N-3][j]
    //            +w/4*(                    0
    //                  + Unew[N-3][j-1]                + Uold[N-3][j+1]
    //                  +Unew[N-4][j])
    //            - w/4*h*h*fun(lam,-2+(j+1)*h,-1+(N-2)*h);
    //        }
    //        // j = M-3;
    //        // k = 0;
    //        Unew[0][M-3]  = (1-w)*Uold[0][M-3]
    //        +w/4*(                Uold[1][M-3]
    //              + Unew[0][M-4]               + Uold[0][M-3]
    //              +    0       )
    //        - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(1)*h);
    //        for (int k=1; k<=N-4; ++k)
    //        {
    //            Unew[k][M-3]  = (1-w)*Uold[k][M-3]
    //            +w/4*(                Uold[k+1][M-3]
    //                  + Unew[k][M-4]              + Uold[k][M-3]
    //                  +Unew[k-1][M-3])
    //            - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(k+1)*h);
    //        }
    //        // k=N-3;
    //        Unew[N-3][M-3] = (1-w)*Uold[N-3][M-3]
    //        +w/4*(                    0
    //              + Unew[N-3] [M-4]               + Uold[N-3][M-3]
    //              +Unew[N-4][M-3])
    //        - w/4*h*h*fun(lam,-2+(M-2)*h,-1+(N-2)*h);
    //
    //        // construct residual matrix res
    //        //j = 0;
    //        //k = 0;
    //        res[0][0]  = -h*h/4*fun(lam,-2+(1)*h,-1+(1)*h) + 0.25*(-4*Unew[0][0]+0+Unew[0][0] + Unew[0][1] + Unew[1][0]);
    //        for (int k=1; k<=N-4; ++k)
    //        {
    //            res[k][0]  = -h*h/4*fun(lam,-2+(1)*h,-1+(k+1)*h) + 0.25*(-4*Unew[k][0]+Unew[k-1][0]+Unew[k][0] + Unew[k][1] + Unew[k+1][0]);
    //        }
    //        //k=N-3;
    //        res[N-3][0]  = -h*h/4*fun(lam,-2+(1)*h,-1+(N-2)*h) + 0.25*(-4*Unew[N-3][0]+Unew[N-4][0]+Unew[N-3][0] + Unew[N-3][1] + 0);
    //        for (int j=1; j<=M-4; ++j)
    //        {
    //            //k=0;
    //            res[0][j]  = -h*h/4*fun(lam,-2+(j+1)*h,-1+(1)*h) + 0.25*(-4*Unew[0][j]+0+Unew[0][j-1] + Unew[0][j+1] + Unew[1][j]);
    //            for (int k=1; k<=N-4; ++k)
    //            {
    //                res[k][j]  = -h*h/4*fun(lam,-2+(j+1)*h,-1+(k+1)*h) + 0.25*(-4*Unew[k][j]+Unew[k-1][j]+Unew[k][j-1] + Unew[k][j+1] + Unew[k+1][j]);
    //            }
    //            // k=N-3;
    //            res[N-3][j]  = -h*h/4*fun(lam,-2+(j+1)*h,-1+(N-2)*h) + 0.25*(-4*Unew[N-3][j]+Unew[N-4][j]+Unew[N-3][j-1] + Unew[N-3][j+1] + 0);
    //        }
    //        // j = M-3;
    //        // k = 0;
    //        res[0][M-3]  = -h*h/4*fun(lam,-2+(M-2)*h,-1+(1)*h) + 0.25*(-4*Unew[0][M-3]+0+Unew[0][M-4] + Unew[0][M-3] + Unew[1][M-3]);
    //        for (int k=1; k<=N-4; ++k)
    //        {
    //            res[k][M-3]  = -h*h/4*fun(lam,-2+(M-2)*h,-1+(k+1)*h) + 0.25*(-4*Unew[k][M-3]+Unew[k-1][M-3]+Unew[k][M-4] + Unew[k][M-3] + Unew[k+1][M-3]);
    //        }
    //        // k=N-3;
    //        res[N-3][M-3]  = -h*h/4*fun(lam,-2+(M-2)*h,-1+(N-2)*h) + 0.25*(-4*Unew[N-3][M-3]+Unew[N-4][M-3]+Unew[N-3][M-4] + Unew[N-3][M-3] + 0);
    //
    //        // compute norm of resisual
    //        //        resnorm =findNorm(&res[0][0],M-2,N-2);
    //        resnorm = findNorm(M-2,N-2,res);
    //        resnormvec[count] = resnorm;
    //
    //
    //
    //        // Uold = Unew - update solution
    //        for (int p=0; p<=M-3; ++p)
    //        {
    //            for (int q=0; q<=N-3; ++q)
    //            {
    //                Uold[q][p] = Unew[q][p];
    //            }
    //        }
    //
    //        count++;
    //    }
    // get ending time
    gettimeofday(&t_end, 0);
    //    double t_elapsed = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)*1e-6; // compute time elapsed in seconds
    //    printf("The total number of iterations is %d \n",count);
    //    printf("The norm of the residual is %5.5e \n",resnorm);
    //    printf("The total elaspsed time is %f\n",t_elapsed);
    //    double t_avg = t_elapsed/(count);
    //    printf("The average time used for a single pass through the grid is %f\n",t_avg);
    //
    //
    //
    // write solution file
    FILE *sol = fopen("sources.out","w");
    fwrite(U,sizeof(double),M*N,sol);
    fclose(sol);
    
    //     //write res matrix
    //    FILE *resmat = fopen("resmat.bin","w");
    //    fwrite(res,sizeof(double),(M-2)*(N-2),resmat);
    //    fclose(resmat);
    
    //    int size = sizeof(resnormvec)/sizeof(resnormvec[0]);
    //    printf("count = %d\n",count);
    //    printf("size of resnormvec = %d\n",size);
    
    // write residual file
    FILE *residual = fopen("residual.out","w");
    fwrite(resmaxvec,sizeof(double),count,residual);
    fclose(residual);
    
    free(resmaxvec);
    free(U);
    
    return 0;
}



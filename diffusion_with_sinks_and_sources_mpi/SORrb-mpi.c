
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

double fun(double lam, double x, double y);
//double findNorm(int M, int N, double res[M][N] );

double fun(double lam, double x, double y)
{
    return  10*lam/sqrt(M_PI)*exp(-lam*lam*((x-1)*(x-1)+y*y))-10*lam/sqrt(M_PI)*exp(-lam*lam*((x+1)*(x+1)+y*y));
}


//double findNorm(int M, int N, double res[N][M])
//{
//    double sum = 0.;
//    for (int m=0; m<M; ++m)
//    {
//        for (int n=0; n<N; ++n)
//        {
//            sum = sum + res[n][m]*res[n][m];
//        }
//    }
//    return sqrt(sum);
//}


int main(int argc, char* argv[])
{
    // initialize the MPI interface
    MPI_Init(&argc, &argv);
    // obtain the rank of current processor
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // obtain the size of the current processor
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // set up parameters
    double w = 1.8;
    int N =128;
    int Ns =N/2+1;
    int M = 2*N-1;
    double tol = 1e-9; // set up tolerance
    double lam = 10;
    double h = 2./(N-1);
    printf("h = %10.10f\n",h);
    int K = 2000;
    int count = 0;
    double* resmaxvec = (double*)malloc(K*sizeof(double)); // maximum residual array
    double resmax = 1000.;
    double val;
    
    // we split the domain half way - each subdomain has (N/2+1)*M grid points
    double* U = (double*)malloc(M*Ns * sizeof(double));
    
    // calculate time - initialize time variables
//    struct timeval t_start;
//    struct timeval t_end;
    
    // initialize U
    for (int j=0; j<=M-1; ++j)
    {
        for (int k=0; k<=Ns-1; ++k)
        {
            U[Ns*j+k] = 0;
        }
    }
    
    // storage for tracking communication info
    MPI_Status status;
    
//    // get starting time
//    gettimeofday(&t_start, 0);
    
    // main iteration loop
    while (count < K && resmax > tol)
    {
        resmax = 0.;
        
        // rank 0: left subdomain
        if (rank == 0)
        {
            for (int j=0; j<=M-1; ++j)
            {
                // receives data from the right subdomain (rank 0 to rank 1)
                MPI_Irecv(&(U[Ns*j+(Ns-1)]), 1, MPI_DOUBLE, 1, j,
                         MPI_COMM_WORLD, &status);
                // sends data to the right subdomain (rank 1 to rank 0)
                MPI_Isend(&(U[Ns-2]), 1, MPI_DOUBLE, 1, j,
                     MPI_COMM_WORLD);
            }
            
        }
        // rank 1: right subdomain
        if (rank == 1)
        {
            // receives data from the left subdomain (rank 1 to rank 0)
            for (int j=0; j<=M-1; ++j)
            {
                // sends data to the left subdomain (rank 1 to rank 0)
                MPI_Isend(&(U[Ns*j+1]), 1, MPI_DOUBLE, rank-1, j, MPI_COMM_WORLD);
                // if data moves right, then rank size-1 receives from the left
                MPI_Irecv(&(U[Ns*j+0]), 1, MPI_DOUBLE, rank-1, j, MPI_COMM_WORLD,
                         &status);
            }
            
        }
        
        
        // now update interior nodes:
        // update red nodes
        // j = 0
        for (int k=2; k<=Ns-2; k = k+2)
        {
                U[Ns*0+k] = (1-w)*U[Ns*0+k]+0.25*w*(U[Ns*0+(k+1)]+ U[Ns*(0)+k]+ U[Ns*(0+1)+k]+U[Ns*0+(k-1)] )- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k= (j+1)%2+1; k<=Ns-2; k = k+2)
            {
                    val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[Ns*j+k]+U[N*j+(k-1)]+U[N*(j-1)+k]+ U[N*(j+1)+k] + U[N*j+(k+1)]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    U[Ns*j+k] = (1-w)*U[Ns*j+k]+0.25*w*(U[Ns*j+(k+1)]+ U[Ns*(j-1)+k]+ U[Ns*(j+1)+k]+U[Ns*j+(k-1)] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
            }
        }
        // j=M-1
        for (int k=2; k<=Ns-2; k=k+2)
        {
                U[Ns*(M-1)+k] = (1-w)*U[Ns*(M-1)+k]+0.25*w*(U[Ns*(M-1)+(k+1)]+ U[Ns*(M-1-1)+k]+ U[Ns*(M-1)+k]+U[Ns*(M-1)+(k-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
        }
        // update black nodes
        // j = 0
        for (int k=1; k<=Ns-2; k=k+2)
        {
                U[Ns*0+k] = (1-w)*U[Ns*0+k]+0.25*w*(U[Ns*0+(k+1)]+ U[Ns*(0)+k]+ U[Ns*(0+1)+k]+U[Ns*0+(k-1)] )- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k=j%2+1; k<=Ns-2; k=k+2)
            {
                    U[Ns*j+k] = (1-w)*U[Ns*j+k]+0.25*w*(U[Ns*j+(k+1)]+ U[Ns*(j-1)+k]+ U[Ns*(j+1)+k]+U[Ns*j+(k-1)] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
                    val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[Ns*j+k]+U[Ns*j+(k-1)]+U[Ns*(j-1)+k]+ U[Ns*(j+1)+k] + U[Ns*j+(k+1)]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
            }
        }
        // j=M-1
        for (int k=1; k<=Ns-2; k=k+2)
        {
                U[Ns*(M-1)+k] = (1-w)*U[Ns*(M-1)+k]+0.25*w*(U[Ns*(M-1)+(k+1)]+ U[Ns*(M-1-1)+k]+ U[Ns*(M-1)+k]+U[Ns*(M-1)+(k-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
        }
        resmaxvec[count] = resmax;
        count++;
    }
    
    // now we gather data by sending all data to rank 0
    if (rank ==0)
    {
           double* finalu = (double*)malloc(N*M * sizeof(double));
        
    }
    
    
    
//    gettimeofday(&t_end, 0);
//        double t_elapsed = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)*1e-6; // compute time elapsed in seconds
    //    printf("The total number of iterations is %d \n",count);
    //    printf("The norm of the residual is %5.5e \n",resnorm);
//        printf("The total elaspsed time is %f\n",t_elapsed);
    //    double t_avg = t_elapsed/(count);
    //    printf("The average time used for a single pass through the grid is %f\n",t_avg);
    //
    //
    //
//    // write solution file
//    FILE *sol = fopen("sources.out","w");
//    fwrite(U,sizeof(double),M*N,sol);
//    fclose(sol);
    
    //     //write res matrix
    //    FILE *resmat = fopen("resmat.bin","w");
    //    fwrite(res,sizeof(double),(M-2)*(N-2),resmat);
    //    fclose(resmat);
    
    //    int size = sizeof(resnormvec)/sizeof(resnormvec[0]);
    //    printf("count = %d\n",count);
    //    printf("size of resnormvec = %d\n",size);
    
//    // write residual file
//    FILE *residual = fopen("residual.out","w");
//    fwrite(resmaxvec,sizeof(double),count,residual);
//    fclose(residual);
    
    free(resmaxvec);
    free(U);
    MPI_Finalize();
    return 0;
}





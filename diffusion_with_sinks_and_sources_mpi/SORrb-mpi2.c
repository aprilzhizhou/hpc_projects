
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
    int Ns =N/2+1; // local N
    int M = 2*N-1;
    double tol = 1e-9; // set up tolerance
    double lam = 10;
    double h = 2./(N-1);
    printf("h = %10.10f\n",h);
    int K = 10;
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
            U[M*k+j] = 0;
        }
    }
    
    
//    // get starting time
//    gettimeofday(&t_start, 0);
    
    // Storage for tracking communicaton info
    MPI_Request sendLeftRequest;
    MPI_Request sendRightRequest;
    MPI_Request recvLeftRequest;
    MPI_Request recvRightRequest;
    
    // main loop
    while (count < K && resmax > tol)
    {
        // rank 1: right subdomain
        if (rank > 1)
        {
            // sends data to the left subdomain (rank 1 to rank 0)
            MPI_Isend(&(U[M*1]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,
                      &sendLeftRequest);
            // receives data from the left subdomain (rank 1 to rank 0)
            MPI_Irecv(&(U[M*0]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,
                      &recvLeftRequest);
        }
        
        // rank 0: left subdomain
        if (rank < size-1)
        {
            // receives data from the right subdomain (rank 0 to rank 1)
            MPI_Isend(&(U[M*(Ns-2)]), M, MPI_DOUBLE, rank+1, count,
                      MPI_COMM_WORLD, &sendRightRequest);
            // sends data to the right subdomain (rank 1 to rank 0)
            MPI_Irecv(&(U[M*(Ns-1)]), M, MPI_DOUBLE, rank+1, count,
                      MPI_COMM_WORLD, &recvRightRequest);
        }
    
        // now update interior nodes:
    resmax = 0;
    // update red nodes
    // j = 0
    for (int k=2; k<=Ns-2; k = k+2)
    {
        U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
    }
    for (int j=1; j<=M-2; ++j)
    {
        for (int k= (j+1)%2+1; k<=Ns-2; k = k+2)
        {
            val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
            if (val > resmax)
            {
                resmax = val;
            }
            U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
        }
    }
    // j=M-1
    for (int k=2; k<=Ns-2; k=k+2)
    {
        U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
    }
    // update black nodes
    // j = 0
    for (int k=1; k<=Ns-2; k=k+2)
    {
        U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fun(lam,-2+0*h,-1+k*h);
    }
    for (int j=1; j<=M-2; ++j)
    {
        for (int k=j%2+1; k<=Ns-2; k=k+2)
        {
            
            U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fun(lam,-2+j*h,-1+k*h);
            val =  -h*h*0.25*fun(lam,-2+j*h,-1+k*h) + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
            if (val > resmax)
            {
                resmax = val;
            }
            //}
        }
    }
    // j=M-1
    for (int k=1; k<=Ns-2; k=k+2)
    {
        U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fun(lam,-2+(M-1)*h,-1+k*h);
        //}
    }
    resmaxvec[count] = resmax;
    count++;
    }


    
    // now we gather data by sending all data to rank 0
    if (rank == 0) {
        // allocate space for the full array
        double* finalu = (double*)malloc(N*M * sizeof(double));
        
        // copy the local data to the full array
        for (int j=0; j<=M-1; ++j )
        {
            for (int k=0; k<=Ns-2; ++k)
            {
                finalu[M*k+j] = U[M*k+j];
            }
        }
//        // Track where the next array data will be inserted in the
//        //      full array
//        int nextj = (Ns-1)*M;
//        int datalen;
//
//        // request the data from each of the other processes,
//        //      appending as we go.
//        for (int r=1; r<size; ++r) {
//            // amount of data in rank=r process excluding internal boundary
//            //      data
//            datalen = N/size + (N%size > r ? 1 : 0);
//            MPI_Recv(&(finalu[nextj]), datalen, MPI_DOUBLE, r, 0,
//                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            nextj += datalen; // update the new insertion point
//        }
//
        // Store the final array in binary format to file finalu.dat
        FILE* fileid = fopen("finalu.bin", "w");
        fwrite(finalu, sizeof(double), N*M, fileid);
        fclose(fileid);
    
        

        // Release the finalu allocation of memory
        free(finalu);

    } else {
        // All ranks other than 0 must send their interior data to rank 0
        // The last rank doesn't have a boundary point,
        //      so it's one bigger than the others.
//        if (rank < size-1)
//            MPI_Send(&(U[M]), (Ns-2)*M, MPI_DOUBLE, 0, 0,
//                     MPI_COMM_WORLD);
//        else
            MPI_Send(&(U[M]), (Ns-1)*M, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
    }
    

//    MPI_Gather(&(U[M]), (Ns-1)*M, MPI_DOUBLE, finalU, (Ns-1)*M, MPI_DOUBLE, 0,  

//
//                     FILE* fileid = fopen("finalu.bin", "w");
//                     fwrite(finalu, sizeof(double), N*M, fileid);
//                     fclose(fileid);

    
    
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
    
    
    
    free(resmaxvec);
    free(U);
    MPI_Finalize();
    return 0;
}






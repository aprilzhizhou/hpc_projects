#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include "mpierror.h"

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
    int N =128; //grid size
    int M = 2*N-1;
    int Ns =N/2+1; // grid size for subdomains
    double tol = 1e-9; // set up tolerance
    double lam = 10;
    double h = 2./(N-1);
    printf("h = %10.10f\n",h);
    int K = 200;
    int count = 0;
    double* resmaxvec = (double*)malloc(K*sizeof(double)); // maximum residual array
    double resmax = 1000.;
    double val;
    double fvec[M][Ns];
    
        if (rank == 0)
        {
            for (int j=0; j<=M-1; ++j)
            {
                for (int k=0; k<=Ns-1; ++k)
                {
                    fvec[j][k] = fun(lam,-2+j*h,-1+k*h);
                }
            }
//            FILE *F = fopen("F.out","w");
//            fwrite(fvec,sizeof(double),Ns*M,F);
        }
        if (rank == 1)
        {
            for (int j=0; j<=M-1; ++j)
            {
                for (int k=0; k<=Ns-1; ++k)
                {
                    fvec[j][k] = fun(lam,-2 + j*h,-1+(Ns-2+k)*h);
                }
            }
            
//            FILE *F = fopen("F.out","w");
//            fwrite(fvec,sizeof(double),Ns*M,F);
        }
    
    
    
    
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
        if (rank == 1)
        {
            // sends data to the left subdomain (rank 1 to rank 0)
            //CheckError( MPI_Send(&data, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD));
            CheckError( MPI_Isend(&(U[M*1]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,
                      &sendLeftRequest));
            // receives data from the left subdomain (rank 1 to rank 0)
            CheckError( MPI_Irecv(&(U[M*0]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,
                      &recvLeftRequest));
        }
        
        // rank 0: left subdomain
        if (rank == 0)
        {
            // receives data from the right subdomain (rank 0 to rank 1)
            CheckError( MPI_Isend(&(U[M*(Ns-2)]), M, MPI_DOUBLE, rank+1, count,
                      MPI_COMM_WORLD, &sendRightRequest));
            // sends data to the right subdomain (rank 1 to rank 0)
            CheckError( MPI_Irecv(&(U[M*(Ns-1)]), M, MPI_DOUBLE, rank+1, count,
                      MPI_COMM_WORLD, &recvRightRequest));
        }
        
        // now update interior nodes: k = 2, ..., Ns-3
        resmax = 0;
        // update red nodes: mod(k+j,2)=0
        // j = 0
        for (int k=2; k<=Ns-3; k = k+2)
        {
            U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fvec[0][k];
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k= j%2+2; k<=Ns-3; k = k+2)
            {
                val =  -h*h*0.25*fvec[j][k] + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
                if (val > resmax)
                {
                    resmax = fabs(val);
                }
                U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fvec[j][k];
            }
        }
        // j=M-1
        for (int k=2; k<=Ns-3; k=k+2)
        {
            U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][k];
        }
        // update black nodes: mod(j+k,2)=1
        // j = 0
        for (int k=3; k<=Ns-3; k=k+2)
        {
            U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fvec[0][k];
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k=(j+1)%2+2; k<=Ns-3; k=k+2)
            {
                
                U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fvec[j][k];
                val =  -h*h*0.25*fvec[j][k] + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
                if (val > resmax)
                {
                    resmax = val;
                }
                //}
            }
        }
        // j=M-1
        for (int k=3; k<=Ns-3; k=k+2)
        {
            U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][k];
        }
        
        // now update points next to the boundary y = +/-1
        if (rank ==0) // left domain (update points k=1)
        {
            // k = 1
            // update red nodes: mod(k+j,2)=0
            for (int j=1; j<=M-2; j=j+2)
            {
                    val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
            }
            // update black nodes: mod(k+j,2)=1
            // j = 0, k=1
                U[M*1+0] = (1-w)*U[M*1+0]+0.25*w*(U[M*(1+1)+0]+ U[M*1+0]+ U[M*1+1]+U[M*(1-1)+0])- 0.25*w*h*h* fvec[0][1];
            for (int j=2; j<=M-2; j=j+2)
            {
                //k=1
                U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
                    val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
            }
            // j=M-1 (even)
            U[M*1+(M-1)] = (1-w)*U[M*1+(M-1)]+0.25*w*(U[M*(1+1)+(M-1)]+ U[M*1+(M-2)]+ U[M*1+(M-1)]+U[M*(1-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][1];
        }
        
        if (rank ==1) // right domain (update points k=(Ns-2) (odd))
        {
            // update red nodes: mod(k+j,2)=0
            for (int j=1; j<=M-2; j=j+2)
            {
                    val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                // k = (Ns-2) odd
                U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
            }
            // update black nodes: mod(k+j,2)=1
            // j = 0
                U[M*(Ns-2)+0] = (1-w)*U[M*(Ns-2)+0]+0.25*w*(U[M*((Ns-2)+1)+0]+ U[M*(Ns-2)+0]+ U[M*(Ns-2)+1]+U[M*((Ns-2)-1)+0])- 0.25*w*h*h* fvec[0][(Ns-2)];
            for (int j=2; j<=M-2; j=j+2)
            {
                    U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
                    val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
            }
            // j=M-1
                U[M*(Ns-2)+(M-1) ] = (1-w)*U[M*(Ns-2)+(M-1)]+0.25*w*(U[M*((Ns-2)+1)+(M-1)]+ U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+(M-1)]+U[M*((Ns-2)-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][(Ns-2)];
        }
        
        // now check if boundary data is ready
        // Indicate whether data is ready {sent left, received left,
        //      sent right, received right}
        int ready[4] = {0, 0, 0, 0};
        // Indicate that the update of the end point has been done after the
        //      data transfer.
        bool done[2] = {false, false};
        
        // There is no data transfer at the ends, so mark those as ready.
        if (rank == 0) {
            ready[0] = 1;
            ready[1] = 1;
        }
        if (rank == size-1) {
            ready[2] = 1;
            ready[3] = 1;
        }
        
        // Keep checking until data is ready and end points updated.
        while (!done[0] || !done[1])
        {
            // Check whether interchange of left data is complete
            if (rank ==1 && !done[0]) {
                if (!ready[0])
                    MPI_Test(&sendLeftRequest, &(ready[0]), MPI_STATUS_IGNORE);
                if (!ready[1])
                    MPI_Test(&recvLeftRequest, &(ready[1]), MPI_STATUS_IGNORE);
            }

            // If the data is exchanged and endpoint hasn't been updated yet,
            //      then update it.
            if (ready[0] && ready[1] && !done[0]) {
                // data on the left has been sent and received, it's safe to
                // update (k = 1)
                // k = 1
                // update red nodes: mod(k+j,2)=0
                for (int j=1; j<=M-2; j=j+2)
                {
                    val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
                }
                // update black nodes: mod(k+j,2)=1
                // j = 0, k=1
                U[M*1+0] = (1-w)*U[M*1+0]+0.25*w*(U[M*(1+1)+0]+ U[M*1+0]+ U[M*1+1]+U[M*(1-1)+0])- 0.25*w*h*h* fvec[0][1];
                for (int j=2; j<=M-2; j=j+2)
                {
                    U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
                    val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                }
                // j=M-1 (even)
                U[M*1+(M-1)] = (1-w)*U[M*1+(M-1)]+0.25*w*(U[M*(1+1)+(M-1)]+ U[M*1+(M-2)]+ U[M*1+(M-1)]+U[M*(1-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][1];
                
                done[0] = true;
            }

            // Check whether interchange of right data is complete
            if (rank ==0 && !done[1]) {
                if (!ready[2])
                    MPI_Test(&sendRightRequest, &(ready[2]), MPI_STATUS_IGNORE);
                if (!ready[3])
                    MPI_Test(&recvRightRequest, &(ready[3]), MPI_STATUS_IGNORE);
            }

            // If the data is exchanged and endpoint hasn't been updated yet,
            //      then update it
            if (ready[2] && ready[3] && !done[1]) {
                // data on the right has been sent and received, it's safe to
                //      update k=Ns-2 (odd)
                
                // update red nodes: mod(k+j,2)=0
                for (int j=1; j<=M-2; j=j+2)
                {
                    val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    // k = (Ns-2) odd
                    U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
                }
                // update black nodes: mod(k+j,2)=1
                // j = 0
                U[M*(Ns-2)+0] = (1-w)*U[M*(Ns-2)+0]+0.25*w*(U[M*((Ns-2)+1)+0]+ U[M*(Ns-2)+0]+ U[M*(Ns-2)+1]+U[M*((Ns-2)-1)+0])- 0.25*w*h*h* fvec[0][(Ns-2)];
                for (int j=2; j<=M-2; j=j+2)
                {
                    U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
                    val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                }
                // j=M-1
                U[M*(Ns-2)+(M-1) ] = (1-w)*U[M*(Ns-2)+(M-1)]+0.25*w*(U[M*((Ns-2)+1)+(M-1)]+ U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+(M-1)]+U[M*((Ns-2)-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][(Ns-2)];
                
                
                done[1] = true;
            }
        }

        
        resmaxvec[count] = resmax;
        count++;
    }
    
    
    // now we gather data by sending all data to rank 0
    double* finalu;
    if (rank == 0) {
        // allocate space for the full array
        finalu = (double*)malloc((N)*M * sizeof(double));
        // copy the local data to the full array
        for (int j=0; j<=M-1; ++j )
        {
            for (int k=0; k<=Ns-2; ++k)
            {
                finalu[M*k+j] = U[M*k+j];
            }
        }
    }
    CheckError( MPI_Gather(&(U[M]),(Ns-1)*M, MPI_DOUBLE, finalu, (Ns-1)*M, MPI_DOUBLE, 0, MPI_COMM_WORLD));

    
    if (rank ==0)
    {
        // Store the final array in binary format to file finalu.dat
        FILE* fileid = fopen("finalu.bin", "w");
        fwrite(finalu, sizeof(double), N*M, fileid);
        fclose(fileid);
        
        // write residual file
        FILE *residual = fopen("residual.out","w");
        fwrite(resmaxvec,sizeof(double),count,residual);
        
        free(finalu);
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
    
    
    
    free(resmaxvec);
    free(U);
    MPI_Finalize();
    return 0;
}







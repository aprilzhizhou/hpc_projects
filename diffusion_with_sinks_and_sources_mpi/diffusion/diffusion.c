#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include "mpierror.h"

double fun(double lam, double x, double y);

double fun(double lam, double x, double y)
{
    double PI = 4.0*atan(1.0);
    return  10*lam/sqrt(PI)*exp(-lam*lam*((x-1)*(x-1)+y*y))-10*lam/sqrt(PI)*exp(-lam*lam*((x+1)*(x+1)+y*y));
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
    
    // read input arguments
    int N = atoi(argv[1]); // grid number
    double w = atof(argv[2]); // relaxation parameter
    double tol = atof(argv[3]); // tolerance
    int K = atoi(argv[4]); // maximum number of iterations
    
    // set up parameters
    int Ns = N/size + (N%size > rank ? 1 : 0) + (rank > 0 ? 1 : 0)
    + (rank < size-1 ? 1 : 0);
    printf("localN = %d\n", Ns);
    int M = 2*N-1;
    double lam = 100;
    double h = 2./(N-1);
    int count = 0;
    double* resmaxvec = (double*)malloc(K*sizeof(double)); // maximum residual array
    double globalresmax = 1000.;
    double resmax;
    double val;
    double fvec[M][Ns];
    
    // define the RHS function array in rank 0 and 1
    if (rank == 0)
    {
        for (int j=0; j<=M-1; ++j)
        {
            for (int k=0; k<=Ns-1; ++k)
            {
                fvec[j][k] = fun(lam,-2+j*h,-1+k*h);
            }
        }
    }
    
    if (rank > 0)
    {
        for (int j=0; j<=M-1; ++j)
        {
            for (int k=0; k<=Ns-1; ++k)
            {
                fvec[j][k] = fun(lam,-2 + j*h,-1+(rank*(N/size)-1+k)*h);
            }
        }
//        // write f
//        FILE* f1 = fopen("fvec1.out","w");
//        fwrite(fvec,sizeof(double),M*Ns,f1);
//        fclose(f1);
    }
    
    
    // we split the domain half way - each subdomain has (N/2+1)*M grid points
    // allocate memory for the solution array U
    double* U = (double*)malloc(M*Ns * sizeof(double));
    
    // initialize U
    for (int j=0; j<=M-1; ++j)
    {
        for (int k=0; k<=Ns-1; ++k)
        {
            U[M*k+j] = 0;
        }
    }
    
    // calculate time - initialize time variables
    double precision = MPI_Wtick();
    double starttime = MPI_Wtime();
    
    
    // Storage for tracking communicaton info
    MPI_Request sendLeftRequest;
    MPI_Request sendRightRequest;
    MPI_Request recvLeftRequest;
    MPI_Request recvRightRequest;
    
    // main loop
    while (count < K && globalresmax > tol)
    {
        // PHASE 1: update red nodes using the black values passed from the other processor from the last iteration
        resmax = 0;
        if (rank < size -1 ) // update using data passed from the right
        {
            // update red values at row k = (Ns-2) odd
            for (int j=1; j<=M-2; j=j+2)
            {
                val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
                U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
            }
        }
        if (rank > 0 ) // update using data passed from the left
        {
            // update red values at row k = 1
            for (int j=1; j<=M-2; j=j+2)
            {
                val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
                U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
            }
        }
        
        
            // PHASE 2: in phase 2, each processor sends and receives red values from the current iteration and use them to compute black values
            // now all the red boundary nodes have been updated and the black boundary nodes have not been updated
            // we send and receive: red boudnary values from the current iteration and black boundary values from the previous iteration
            
            //  exchange end data to the left
            if (rank > 0)
            {
                // send and receive boundary data from rank 1 to rank 0
                CheckError( MPI_Isend(&(U[M*1]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,&sendLeftRequest));
                CheckError( MPI_Irecv(&(U[M*0]), M, MPI_DOUBLE, rank-1, count, MPI_COMM_WORLD,&recvLeftRequest));
            }
            
            // exchange end data to the right
            if (rank < size-1)
            {
                // send and receive boundary data from rank 0 to rank 1
                CheckError( MPI_Isend(&(U[M*(Ns-2)]), M, MPI_DOUBLE, rank+1, count,MPI_COMM_WORLD, &sendRightRequest));
                CheckError( MPI_Irecv(&(U[M*(Ns-1)]), M, MPI_DOUBLE, rank+1, count,MPI_COMM_WORLD, &recvRightRequest));
            }
        
        // update interior nodes: row k = 2, ..., Ns-3
        // first update red nodes where mod(k+j,2)=0
        // j = 0
        for (int k=2; k<=Ns-3; k = k+2)
        {
            val =  -h*h*0.25*fvec[0][k] + 0.25*(-4*U[M*k+0]+U[M*(k-1)+0]+U[M*k+(0)]+ U[M*k+(1)]+ U[M*(k+1)+0]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
            U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fvec[0][k];
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k= j%2+2; k<=Ns-3; k = k+2)
            {
                val =  -h*h*0.25*fvec[j][k] + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
                U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fvec[j][k];
            }
        }
        // j=M-1
        for (int k=2; k<=Ns-3; k=k+2)
        {
            val =  -h*h*0.25*fvec[M-1][k] + 0.25*(-4*U[M*k+M-1]+U[M*(k-1)+M-1]+U[M*k+(M-1-1)]+ U[M*k+(M-1)]+ U[M*(k+1)+M-1]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
            U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][k];
        }
        // update black nodes: mod(j+k,2)=1
        // j = 0
        for (int k=3; k<=Ns-3; k=k+2)
        {
            U[M*k+0] = (1-w)*U[M*k+0]+0.25*w*(U[M*(k+1)+0]+ U[M*k+0]+ U[M*k+1]+U[M*(k-1)+0])- 0.25*w*h*h* fvec[0][k];
            val =  -h*h*0.25*fvec[0][k] + 0.25*(-4*U[M*k+0]+U[M*(k-1)]+U[M*k]+ U[M*k+1]+ U[M*(k+1)]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
            
        }
        for (int j=1; j<=M-2; ++j)
        {
            for (int k=(j+1)%2+2; k<=Ns-3; k=k+2)
            {
                
                U[M*k+j]= (1-w)*U[M*k+j]+0.25*w*(U[M*(k+1)+j]+ U[M*k+(j-1)]+ U[M*k+(j+1)]+U[M*(k-1)+j] )- 0.25*w*h*h* fvec[j][k];
                val =  -h*h*0.25*fvec[j][k] + 0.25*(-4*U[M*k+j]+U[M*(k-1)+j]+U[M*k+(j-1)]+ U[M*k+(j+1)]+ U[M*(k+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
            }
        }
        // j=M-1
        for (int k=3; k<=Ns-3; k=k+2)
        {
            U[M*k+(M-1)] = (1-w)*U[M*k+(M-1)]+0.25*w*(U[M*(k+1)+(M-1)]+ U[M*k+(M-2)]+ U[M*k+(M-1)]+U[M*(k-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][k];
            val =  -h*h*0.25*fvec[M-1][k] + 0.25*(-4*U[M*k+M-1]+U[M*(k-1)+M-1]+U[M*k+(M-2)]+ U[M*k+(M-1)]+ U[M*(k+1)+M-1]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
        }
        
        // now update boundary nodes that do not depend on data from the other processor
        if (rank ==0) // left domain (update row k=1)
        {
            // k = 1
            // update red nodes: mod(k+j,2)=0
            for (int j=1; j<=M-2; j=j+2)
            {
                val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
                U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
            }
            // update black nodes: mod(k+j,2)=1
            // j = 0, k=1
            U[M*1+0] = (1-w)*U[M*1+0]+0.25*w*(U[M*(1+1)+0]+ U[M*1+0]+ U[M*1+1]+U[M*(1-1)+0])- 0.25*w*h*h* fvec[0][1];
            val =  -h*h*0.25*fvec[0][1] + 0.25*(-4*U[M*1]+U[M*(1-1)]+U[M*1]+ U[M*1+(1)]+ U[M*(1+1)]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
            for (int j=2; j<=M-2; j=j+2)
            {
                //k=1
                U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
                val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
            }
            // j=M-1 (even)
            U[M*1+(M-1)] = (1-w)*U[M*1+(M-1)]+0.25*w*(U[M*(1+1)+(M-1)]+ U[M*1+(M-2)]+ U[M*1+(M-1)]+U[M*(1-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][1];
            val =  -h*h*0.25*fvec[M-1][1] + 0.25*(-4*U[M*1+M-1]+U[M*(1-1)+M-1]+U[M*1+(M-2)]+ U[M*1+(M-1)]+ U[M*(1+1)+M-1]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
        }
        
        if (rank == size-1) // right domain (update row k=(Ns-2) (odd))
        {
            // update red nodes: mod(k+j,2)=0
            for (int j=1; j<=M-2; j=j+2)
            {
                val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                val = fabs(val);
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
            val =  -h*h*0.25*fvec[0][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)]+U[M*((Ns-2)-1)]+U[M*(Ns-2)]+ U[M*(Ns-2)+(1)]+ U[M*((Ns-2)+1)]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
            for (int j=2; j<=M-2; j=j+2)
            {
                U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
                val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                val = fabs(val);
                if (val > resmax)
                {
                    resmax = val;
                }
            }
            // j=M-1
            U[M*(Ns-2)+(M-1) ] = (1-w)*U[M*(Ns-2)+(M-1)]+0.25*w*(U[M*((Ns-2)+1)+(M-1)]+ U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+(M-1)]+U[M*((Ns-2)-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][(Ns-2)];
            val =  -h*h*0.25*fvec[M-1][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+M-1]+U[M*((Ns-2)-1)+M-1]+U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+M-1]+ U[M*((Ns-2)+1)+M-1]);
            val = fabs(val);
            if (val > resmax)
            {
                resmax = val;
            }
        }
        
        
            // now check if the boundary data has been sent and received, and update black nodes using the updated red values receivd
            // Indicate whether data is ready {sent left, received left,
            //      sent right, received right}
            int ready[4] = {0, 0, 0, 0};
            // Indicate that the update of the end point has been done after the
            //      data transfer.
            bool done[2] = {false, false};
            
            // There is no data transfer at the ends, so mark those as ready.
            if (rank == 0) { // rank = 0
                ready[0] = 1;
                ready[1] = 1;
            }
            if (rank == size-1) { // rank = 1
                ready[2] = 1;
                ready[3] = 1;
            }
            
            // Keep checking until data is ready and end points updated.
            while (!done[0] || !done[1])
            {
                // Check whether interchange of left data is complete
                if (rank > 0 && !done[0]) {
                    if (!ready[0])
                        MPI_Test(&sendLeftRequest, &(ready[0]), MPI_STATUS_IGNORE);
                    if (!ready[1])
                        MPI_Test(&recvLeftRequest, &(ready[1]), MPI_STATUS_IGNORE);
                }
                // If the data is exchanged and endpoint hasn't been updated yet,
                //      then update it.
                // now update boudnary red values using the previous black values (from the other processor)
                if (ready[0] && ready[1] && !done[0]) {
                    // data on the left has been sent and received, it's safe to
                    // update nodes on row k = 1
                    // update black nodes where mod(k+j,2)=1
                    // j = 0, k=1
                    U[M*1+0] = (1-w)*U[M*1+0]+0.25*w*(U[M*(1+1)+0]+ U[M*1+0]+ U[M*1+1]+U[M*(1-1)+0])- 0.25*w*h*h* fvec[0][1];
                    val =  -h*h*0.25*fvec[0][1] + 0.25*(-4*U[M*1]+U[M*(1-1)]+U[M*1]+ U[M*1+(1)]+ U[M*(1+1)]);
                    val = fabs(val);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    for (int j=2; j<=M-2; j=j+2)
                    {
                        U[M*1+j]= (1-w)*U[M*1+j]+0.25*w*(U[M*(1+1)+j]+ U[M*1+(j-1)]+ U[M*1+(j+1)]+U[M*(1-1)+j] )- 0.25*w*h*h* fvec[j][1];
                       val =  -h*h*0.25*fvec[j][1] + 0.25*(-4*U[M*1+j]+U[M*(1-1)+j]+U[M*1+(j-1)]+ U[M*1+(j+1)]+ U[M*(1+1)+j]);
                       val = fabs(val);
                       if (val > resmax)
                       {
                           resmax = val;
                       }
                    }
                    // j=M-1 (even)
                    U[M*1+(M-1)] = (1-w)*U[M*1+(M-1)]+0.25*w*(U[M*(1+1)+(M-1)]+ U[M*1+(M-2)]+ U[M*1+(M-1)]+U[M*(1-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][1];
                    val =  -h*h*0.25*fvec[M-1][1] + 0.25*(-4*U[M*1+M-1]+U[M*(1-1)+M-1]+U[M*1+(M-2)]+ U[M*1+(M-1)]+ U[M*(1+1)+M-1]);
                    val = fabs(val);
                    if (val > resmax)
                    {
                        resmax = val;
                    }
                    
                        done[0] = true;
                }
                    
                    // Check whether interchange of right data is complete
                    if (rank < size-1 && !done[1]) {
                        if (!ready[2])
                            MPI_Test(&sendRightRequest, &(ready[2]), MPI_STATUS_IGNORE);
                        if (!ready[3])
                            MPI_Test(&recvRightRequest, &(ready[3]), MPI_STATUS_IGNORE);
                    }
                    
                    // If the data is exchanged and endpoint hasn't been updated yet,
                    //      then update it
                    // now update boudnary red values using the previous black values (from the other processor)
                    
                    if (ready[2] && ready[3] && !done[1]) {
                        // data on the right has been sent and received, it's safe to
                        // update nodes on row k=Ns-2 (odd)
                        // update black nodes where mod(k+j,2)=1
                        // j = 0
                        U[M*(Ns-2)+0] = (1-w)*U[M*(Ns-2)+0]+0.25*w*(U[M*((Ns-2)+1)+0]+ U[M*(Ns-2)+0]+ U[M*(Ns-2)+1]+U[M*((Ns-2)-1)+0])- 0.25*w*h*h* fvec[0][(Ns-2)];
                        val =  -h*h*0.25*fvec[0][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)]+U[M*((Ns-2)-1)]+U[M*(Ns-2)]+ U[M*(Ns-2)+(1)]+ U[M*((Ns-2)+1)]);
                        val = fabs(val);
                        if (val > resmax)
                        {
                            resmax = val;
                        }
                        for (int j=2; j<=M-2; j=j+2)
                        {
                            U[M*(Ns-2)+j]= (1-w)*U[M*(Ns-2)+j]+0.25*w*(U[M*((Ns-2)+1)+j]+ U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+U[M*((Ns-2)-1)+j] )- 0.25*w*h*h* fvec[j][(Ns-2)];
                            val =  -h*h*0.25*fvec[j][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+j]+U[M*((Ns-2)-1)+j]+U[M*(Ns-2)+(j-1)]+ U[M*(Ns-2)+(j+1)]+ U[M*((Ns-2)+1)+j]);
                            val = fabs(val);
                            if (val > resmax)
                            {
                                resmax = val;
                            }
                        }
                        // j=M-1
                        U[M*(Ns-2)+(M-1) ] = (1-w)*U[M*(Ns-2)+(M-1)]+0.25*w*(U[M*((Ns-2)+1)+(M-1)]+ U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+(M-1)]+U[M*((Ns-2)-1)+(M-1)] )- 0.25*w*h*h* fvec[M-1][(Ns-2)];
                        val =  -h*h*0.25*fvec[M-1][(Ns-2)] + 0.25*(-4*U[M*(Ns-2)+M-1]+U[M*((Ns-2)-1)+M-1]+U[M*(Ns-2)+(M-2)]+ U[M*(Ns-2)+(M-1)]+ U[M*((Ns-2)+1)+M-1]);
                        val = fabs(val);
                        if (val > resmax)
                        {
                            resmax = val;
                        }
                        
                        done[1] = true;
                    }
                }

            
        // obtain maximum residual from both processors
        CheckError( MPI_Allreduce(&resmax, &globalresmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD));
        
        resmaxvec[count] = globalresmax;
        count++;
    }
    
    // now we gather data by sending all data to rank 0
    double* finalu;
    if (rank == 0) {
        // allocate space for the full array
        finalu = (double*)malloc(N*M * sizeof(double));
        // copy the local data to the full array
        for (int j=0; j<=M-1; ++j )
        {
            for (int k=0; k<=Ns-2; ++k)
            {
                finalu[M*k+j] = U[M*k+j];
            }
        }
    }
    
    
    // get finishing time
    double time_elapsed = MPI_Wtime() - starttime;
    printf("Execution time = %le seconds, with precision %le seconds\n",time_elapsed, precision);
    
    // collect final solution
    CheckError( MPI_Gather(&(U[M]),(N/size)*M, MPI_DOUBLE, finalu, (N/size)*M, MPI_DOUBLE, 0, MPI_COMM_WORLD));
    
    
    if (rank ==0)
    {
        printf("The total number of iterations is %d \n",count);
        printf("The the maximum residual is %5.5e \n",globalresmax);
        
        // Store the final array in binary format to file finalu.dat
        FILE* fileid = fopen("sources.out", "w");
        fwrite(finalu, sizeof(double), N*M, fileid);
        fclose(fileid);
        
//        // write residual file
//        FILE *residual = fopen("residual.out","w");
//        fwrite(resmaxvec,sizeof(double),count,residual);
        free(finalu);
    }
    
    
    free(resmaxvec);
    free(U);
    MPI_Finalize();
    return 0;
}







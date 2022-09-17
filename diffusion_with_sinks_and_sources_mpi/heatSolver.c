#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>  // This is the header file for MPI functions

/*
 int main(int argc, char* argv[])
 
 Heat flow solver.  All data is initially zero, boundary conditions are
 provided as input
 
 Inputs: argc should be 2
 argv[1]: number of grid points
 User prompted for left and right boundary conditions
 
 Outputs: Prints out the final values.
 
 */

int main(int argc, char* argv[])
{
    // First thing is to initialize the MPI interface.
    // Some arguments can be passed through to the MPI interface,
    // so the argument list is sent by sending the argument list
    MPI_Init(&argc, &argv);
    
    int N = atoi(argv[1]); // Get the number of grid points
    
    // Determine the rank of the current process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Determine the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // allocate memory.  Each core should have roughly N/size + 2 grid
    // points, but if N is not evenly divisible by size, give one extra
    // grid point to the low ranks to make up the difference.  Also, the
    // rank=0 and rank=size-1 processes do not have to store boundary
    // data for a neighbor.
    int localN = N/size + (N%size > rank ? 1 : 0) + (rank > 0 ? 1 : 0)
    + (rank < size-1 ? 1 : 0);
    double* u[2];
    // need two copies to do updates
    u[0] = (double*)malloc(localN * sizeof(double));
    u[1] = (double*)malloc(localN * sizeof(double));
    
    // Initialize the data
    for (int i=0; i<localN; ++i)
        u[0][i] = 0.;
    
    if (rank == 0) {
        // Read the boundary conditions from the input list.
        // Get the left boundary condition.
        u[0][0] = atof(argv[2]);
        u[1][0] = u[0][0];
    }
    
    if (rank == size-1) {
        // Get the right boundary condition
        u[0][localN-1] = atof(argv[3]);
        u[1][localN-1] = u[0][localN-1];
    }
    
    // CFL condition: dt < 0.5 dx^2
    double dx = 1./(N-1);
    double dx2 = dx*dx;
    double dt = 0.25 * dx2;
    
    // Storage for tracking communicaton info
    MPI_Request sendLeftRequest;
    MPI_Request sendRightRequest;
    MPI_Request recvLeftRequest;
    MPI_Request recvRightRequest;
    
    // main loop: terminal time is T=1
    int i, newi, oldi;
    for (i=0; i*dt < 1.0; ++i) {
        newi = (i+1)%2;
        oldi = i%2;
        // Exchange end data to the left,
        //      the tag will correspond to the step #.
        if (rank == 1) {
            MPI_Isend(&(u[oldi][1]), 1, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD,
                      &sendLeftRequest);
            MPI_Irecv(&(u[oldi][0]), 1, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD,
                      &recvLeftRequest);
        }
        // Exchange end data to the right,
        //      the tag will correspond to the step #
        if (rank == 0) {
            MPI_Isend(&(u[oldi][localN-2]), 1, MPI_DOUBLE, rank+1, i,
                      MPI_COMM_WORLD, &sendRightRequest);
            MPI_Irecv(&(u[oldi][localN-1]), 1, MPI_DOUBLE, rank+1, i,
                      MPI_COMM_WORLD, &recvRightRequest);
        }
        
        // Now do update in the interior where it doesn't matter if we have
        //      the current data from neighboring processes
        for (int j=2; j<localN-2; ++j)
            u[newi][j] = u[oldi][j]+dt*(u[oldi][j+1] - 2*u[oldi][j]
                                        + u[oldi][j-1])/dx2;
        
        // Can update points next to boundary
        if (rank == 0)
            u[newi][1] = u[oldi][1]+dt*(u[oldi][2] - 2*u[oldi][1]
                                        + u[oldi][0])/dx2;
        if (rank == size-1)
            u[newi][localN-2] = u[oldi][localN-2]+dt*(u[oldi][localN-1]
                                                      - 2*u[oldi][localN-2]  + u[oldi][localN-3])/dx2;
        
        // Now check to see boundary data is ready.
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
            if (rank > 0 && !done[0]) {
                if (!ready[0])
                    MPI_Test(&sendLeftRequest, &(ready[0]), MPI_STATUS_IGNORE);
                if (!ready[1])
                    MPI_Test(&recvLeftRequest, &(ready[1]), MPI_STATUS_IGNORE);
            }
            
            // If the data is exchanged and endpoint hasn't been updated yet,
            //      then update it.
            if (ready[0] && ready[1] && !done[0]) {
                // data on the left has been sent and received, it's safe to
                //      update now
                u[newi][1] = u[oldi][1] + dt*(u[oldi][2] - 2*u[oldi][1]
                                              + u[oldi][0])/dx2;
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
            if (ready[2] && ready[3] && !done[1]) {
                // data on the right has been sent and received, it's safe to
                //      update now
                u[newi][localN-2] = u[oldi][localN-2] + dt*(u[oldi][localN-1]
                                                            - 2*u[oldi][localN-2] + u[oldi][localN-3])/dx2;
                done[1] = true;
            }
        }
    }
    
    // We have to gather all the data from the various
    //      processes so it can be printed.
    // We'll send all the data to rank==0.
    // We'll see a better way to do this later.
    if (rank == 0) {
        // allocate space for the full array
        double* finalu = (double*)malloc(N * sizeof(double));
        
        // copy the local data to the full array
        for (int j=0; j<localN-1; ++j)
            finalu[j] = u[newi][j];
        
        // Track where the next array data will be inserted in the
        //      full array
        int nextj = localN-1;
        int datalen;
        
        // request the data from each of the other processes,
        //      appending as we go.
        for (int r=1; r<size; ++r) {
            // amount of data in rank=r process excluding internal boundary
            //      data
            datalen = N/size + (N%size > r ? 1 : 0);
            MPI_Recv(&(finalu[nextj]), datalen, MPI_DOUBLE, r, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            nextj += datalen; // update the new insertion point
        }
        
        // Store the final array in binary format to file finalu.dat
        FILE* fileid = fopen("finalu.bin", "w");
        fwrite(finalu, sizeof(double), N, fileid);
        fclose(fileid);
        
        // Release the finalu allocation of memory
        free(finalu);
        
    } else {
        // All ranks other than 0 must send their interior data to rank 0
        // The last rank doesn't have a boundary point,
        //      so it's one bigger than the others.
        if (rank < size-1)
            MPI_Send(&(u[oldi][1]), localN-2, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
        else
            MPI_Send(&(u[oldi][1]), localN-1, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
    }
    
    free(u[0]);
    free(u[1]);
    
    // Must shut down the MPI system when you're done.
    // If you don't, there can be lots of junk left behind.
    MPI_Finalize();
    
    return 0;
}

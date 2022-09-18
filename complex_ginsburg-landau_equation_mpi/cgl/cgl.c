#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <complex.h>
#include <sys/time.h>
#include <time.h>


int main(int argc, char* argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Initialize FFTW for MPI
    fftw_mpi_init();
    
    #ifndef M_PI // define pi
        const double M_PI = 4.0*atan(1.0);
    #endif
    
    // read input arguments
    const ptrdiff_t N = atoi(argv[1]); // grid number
    double c1 = atof(argv[2]); // c1
    double c3 = atof(argv[3]); // c3
    int M = atoi(argv[4]); // number of time steps
    
    
    // Get dimensions of the domain
//    const ptrdiff_t N = 128;
    ptrdiff_t localN, local0;
    
//    int M;
    double L, T, dx, dy, dt;
    dx = 2.*M_PI/N; // space grid size
    dy = 2.*M_PI/N;
    L = 128*M_PI;
    T= 10000.; // terminal time
    dt = T/M; // time step size

    
    // Determine the amount of local memory required
    ptrdiff_t alloc_local = fftw_mpi_local_size_2d(N, N, MPI_COMM_WORLD,&localN, &local0);
    double localNprint = 1.*localN;
    printf("localN = %f\n",localNprint);
    
    // initialize input and output variables and Allocate the local memory
    fftw_complex* A = fftw_alloc_complex(alloc_local);
    fftw_complex* Af = fftw_alloc_complex(alloc_local);
    fftw_complex* Anewf = fftw_alloc_complex(alloc_local);
    fftw_complex* lapf = fftw_alloc_complex(alloc_local);
    fftw_complex* ag = fftw_alloc_complex(alloc_local);
    fftw_complex* agf = fftw_alloc_complex(alloc_local);
    fftw_complex* nl = fftw_alloc_complex(alloc_local);
    fftw_complex* nlf = fftw_alloc_complex(alloc_local);
    fftw_complex* Astore = fftw_alloc_complex(alloc_local);
    
    
    fftw_plan pf, pb_nl, pf_nl, pb, pb_store; // initialize transform plans
    // set up transform plans
    pf = fftw_mpi_plan_dft_2d(N,N, A, Af,MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); // plan to fourier transform initial condition
    pb_nl = fftw_mpi_plan_dft_2d(N,N, agf, ag,MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);// inverse transform for the nonlinear term
    pf_nl = fftw_mpi_plan_dft_2d(N,N, nl, nlf,MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); // Fourier transform for the nonlinear term
    pb = fftw_mpi_plan_dft_2d(N,N, Anewf, A, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); // plan to inverse transform (the final time solution)
    pb_store = fftw_mpi_plan_dft_2d(N,N, Af, Astore, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); // plan to inverse transform (the final time solution)
    
    
    // create random initial values (locally)
    long int seed = (long int)time(NULL); // seed for random initial values
    
    if (argc >5 )
        seed = atol(argv[1]);
    
    printf("Starting seed = %ld\n",seed);
    srand48(seed);
    double rand1;
    double rand2;
    double x;
    double y;
    for (int i=0; i<localN; ++i) {
        for (int j=0; j<N; ++j) {
            rand1 = -1.5+3.*drand48();
            rand2 = -1.5+3.*drand48();
            A[i*N+j][0] = rand1;
            A[i*N+j][1] = rand2;
        }
    }
    
    // transform initial condition from physical space to the Fourier space:
    // A -> Af
    fftw_execute(pf); // fourier transform

    // initialize variable Astorefinal
    int rank;
    double* Astorefinal = NULL;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        Astorefinal = (double*)malloc(2*N*N*sizeof(double));
    }
    
    
     FILE* fileid1 = fopen("CGL.out", "w");
    // measure time
    double precision = MPI_Wtick();
    double starttime = MPI_Wtime();
    
    // step forward in time in Fourier space using RK4
    for (int m =0; m <= M; ++m) {
        // 1. RK4 STAGE 1
        // 1.1. assign Af to argument agf, and compute laplacian A in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Af[i*N+j][0];
                agf[i*N+j][1] = Af[i*N+j][1];
                if (globali < N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali + j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali + j*j)*agf[i*N+j][1];
                }
                if (globali < N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][1];
                }
                if (globali >= N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][1];
                }
                if (globali >= N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][1];
                }
            }
        }
        // 1.2. inverse tranfrom argument agf to physical space: agf -> ag
        fftw_execute(pb_nl);
        // rescale
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                ag[i*N+j][0] = ag[i*N+j][0]/(N*N);
                ag[i*N+j][1] = ag[i*N+j][1]/(N*N);
            }
        }
        // 1.3. compute nonlinear term |A|^2 A in physical space
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                nl[i*N+j][0] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][0];
                nl[i*N+j][1] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][1];
            }
        }
        // 1.4. transform nonlinear term nl to Fourier space: nl -> nlf
        fftw_execute(pf_nl);
        // 1.5. update Anewf at first stage
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Anewf[i*N+j][0] = Af[i*N+j][0] + dt/4.*( agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]));
                Anewf[i*N+j][1] = Af[i*N+j][1] + dt/4.*( agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]));
            }
        }

        // ------------------------------------------------------------------------------------------------------------------
        // 2. RK4 STAGE 2
        // 2.1. assign Anewf to argument arf, compute laplacian Anew in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Anewf[i*N+j][0];
                agf[i*N+j][1] = Anewf[i*N+j][1];
                if (globali < N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali + j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali + j*j)*agf[i*N+j][1];
                }
                if (globali < N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][1];
                }
                if (globali >= N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][1];
                }
                if (globali >= N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][1];
                }
            }
        }
        
        // 2.2. inverse tranfrom argument agf to physical space: agf -> ag
        fftw_execute(pb_nl);
        // rescale
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                ag[i*N+j][0] = ag[i*N+j][0]/(N*N);
                ag[i*N+j][1] = ag[i*N+j][1]/(N*N);
            }
        }
        // 2.3. compute nonlinear term |A+k1*dt/2|^2 (A+k1*dt/2) in physical space
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                nl[i*N+j][0] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][0];
                nl[i*N+j][1] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][1];
            }
        }
        // 2.4. transform nonlinear term nl to Fourier space: nl -> nlf
        fftw_execute(pf_nl);
        // 2.5. update Anewf at second stage
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Anewf[i*N+j][0] = Af[i*N+j][0] + dt/3.*( agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]));
                Anewf[i*N+j][1] = Af[i*N+j][1] +dt/3.*( agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]));
            }
        }
        // ------------------------------------------------------------------------------------------------------------------
        // 3. RK4 STAGE 3
        // 3.1. assign Anewf to argument arf, compute laplacian Anew in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Anewf[i*N+j][0];
                agf[i*N+j][1] = Anewf[i*N+j][1];
                if (globali < N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali + j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali + j*j)*agf[i*N+j][1];
                }
                if (globali < N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][1];
                }
                if (globali >= N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][1];
                }
                if (globali >= N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][1];
                }
            }
        }
        // 3.2. inverse tranfrom argument agf to physical space: agf -> ag
        fftw_execute(pb_nl);
        // rescale
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                ag[i*N+j][0] = ag[i*N+j][0]/(N*N);
                ag[i*N+j][1] = ag[i*N+j][1]/(N*N);
            }
        }
        // 3.3. compute nonlinear term |A+k2*dt/2|^2 (A+k2*dt/2) in physical space
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                nl[i*N+j][0] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][0];
                nl[i*N+j][1] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][1];
            }
        }
        // 3.4. transform nonlinear term nl to Fourier space: nl -> nlf
        fftw_execute(pf_nl);
        // 3.5. update Anewf at third stage
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Anewf[i*N+j][0] = Af[i*N+j][0] + dt/2.*( agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]));
                Anewf[i*N+j][1] = Af[i*N+j][1] +dt/2.*( agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]));
            }
        }
        
        // ------------------------------------------------------------------------------------------------------------------
        // 4. RK4 STAGE 4
        // 4.1. assign Anewf to argument arf, compute laplacian Anew in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Anewf[i*N+j][0];
                agf[i*N+j][1] = Anewf[i*N+j][1];
                if (globali < N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali + j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali + j*j)*agf[i*N+j][1];
                }
                if (globali < N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -(globali*globali+(N-j)*(N-j))*agf[i*N+j][1];
                }
                if (globali >= N/2 && j < N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+j*j)*agf[i*N+j][1];
                }
                if (globali >= N/2 && j >= N/2)
                {
                    lapf[i*N+j][0] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][0];
                    lapf[i*N+j][1] = -((N-globali)*(N-globali)+(N-j)*(N-j))*agf[i*N+j][1];
                }
            }
        }
        // 4.2. inverse tranfrom argument agf to physical space: agf -> ag
        fftw_execute(pb_nl);
        // rescale
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                ag[i*N+j][0] = ag[i*N+j][0]/(N*N);
                ag[i*N+j][1] = ag[i*N+j][1]/(N*N);
            }
        }
        // 4.3. compute nonlinear term |A+k3*dt|^2 (A+k3*dt) in physical space
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                nl[i*N+j][0] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][0];
                nl[i*N+j][1] = (ag[i*N+j][0]*ag[i*N+j][0] + ag[i*N+j][1]*ag[i*N+j][1])*ag[i*N+j][1];
            }
        }
        // 4.4. transform nonlinear term nl to Fourier space: nl -> nlf
        fftw_execute(pf_nl);
        // 4.5. update Anewf at fourth stage
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Anewf[i*N+j][0] = Af[i*N+j][0] + dt*( agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]));
                Anewf[i*N+j][1] = Af[i*N+j][1] +dt*( agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]));
            }
        }
        
        // ------------------------------------------------------------------------------------------------------------------
        
        // 5. replace Af with Anewf
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Af[i*N+j][0] = Anewf[i*N+j][0];
                Af[i*N+j][1] = Anewf[i*N+j][1];
            }
        }
        
        // ------------------------------------------------------------------------------------------------------------------
        // 6. write output file CGL.out for m = 1000*k, k = 0,1,2,..,10

        if ( (m%10000)  == 0 ){
            fftw_execute(pb_store);
            // rescale
            for (int i=0; i<localN; ++i) {
                for (int j=0; j<N; ++j) {
                    Astore[i*N+j][0] = Astore[i*N+j][0]/(N*N);
                    Astore[i*N+j][1] = Astore[i*N+j][1]/(N*N);
                }
            }
            printf("t = %f\n",m*dt);
            MPI_Gather(Astore, 2*localN*N, MPI_DOUBLE, Astorefinal, 2*localN*N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // write sollution
            fwrite(Astorefinal, sizeof(fftw_complex), N*N, fileid1);
        }
        
        
    } // end of main temporal loop
    fclose(fileid1);

    // get ending time
    double time_elapsed = MPI_Wtime() - starttime;
    printf("Execution time = %le seconds, with precision %le seconds\n",time_elapsed, precision);
    
    //we have reached the terminal time in the Fourier space, we inverse transform the coefficients Af back to the physical space
    fftw_execute(pb); // inverse transform: Anewf -> A
    // rescale
    for (int i=0; i<localN; ++i) {
        for (int j=0; j<N; ++j) {
            A[i*N+j][0] = A[i*N+j][0]/(N*N);
            A[i*N+j][1] = A[i*N+j][1]/(N*N);
        }
    }
    
    // Set up rank 0 to receive the full array upon completion
    double* fullA = NULL;
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fullA = (double*)malloc(2*N*N*sizeof(double));
    }
    
    MPI_Gather(A, 2*localN*N, MPI_DOUBLE, fullA, 2*localN*N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // write final solution
//    if (rank == 0) {
//        // write sollution A to file
//        FILE* fileid2 = fopen("finalsol.bin", "w");
//        fwrite(fullA, sizeof(fftw_complex), N*N*2, fileid2);
//        fclose(fileid2);
//    }

    
    fftw_destroy_plan(pb);
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb_nl);
    fftw_destroy_plan(pf_nl);
    fftw_destroy_plan(pb_store);
    
    if (rank == 0) {
        free(fullA);
        free(Astorefinal);
    }
    
    fftw_free(A);
    fftw_free(Af);
    fftw_free(Anewf);
    fftw_free(Astore);
    fftw_free(ag);
    fftw_free(agf);
    fftw_free(nl);
    fftw_free(nlf);
    MPI_Finalize();
    return 0;
}

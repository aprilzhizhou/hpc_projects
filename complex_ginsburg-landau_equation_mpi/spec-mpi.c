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
    
    // Get dimensions of the domain
    const ptrdiff_t N = 64;
    ptrdiff_t localN, local0;
    
    int M;
    double c1, c3, L, T, dx, dy, dt;
    M = 10000; // number of time steps
    dx = 2.*M_PI/N; // space grid size
    dy = 2.*M_PI/N;
    c1 = 1.5; // assign parameter values
    c3 = 0.25;
    L = 64.*M_PI;
    T= 1000.; // terminal time
    dt = T/M; // time step size
    printf("dt = %f\n", dt);

    
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
    fftw_complex* k1 = fftw_alloc_complex(alloc_local);
    fftw_complex* k2 = fftw_alloc_complex(alloc_local);
    fftw_complex* k3 = fftw_alloc_complex(alloc_local);
    fftw_complex* k4 = fftw_alloc_complex(alloc_local);
    
    fftw_plan pf, pb_nl, pf_nl, pb; // initialize transform plans
    // set up transform plans
    pf = fftw_mpi_plan_dft_2d(N,N, A, Af,MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); // plan to fourier transform initial condition
    pb_nl = fftw_mpi_plan_dft_2d(N,N, agf, ag,MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);// inverse transform for the nonlinear term
    pf_nl = fftw_mpi_plan_dft_2d(N,N, nl, nlf,MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); // Fourier transform for the nonlinear term
    pb = fftw_mpi_plan_dft_2d(N,N, Anewf, A, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); // plan to inverse transform (the final time solution)
    
    
    // create random initial values (locally)
    long int seed = (long int)time(NULL); // seed for random initial values
    
    if (argc >1 )
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
    
    
    // calculate time - initialize time variables
    struct timeval t_start;
    struct timeval t_end;
    // get starting time
    gettimeofday(&t_start, 0);
    
    // step forward in time in Fourier space using RK4
    for (int m =0; m <= M; ++m) {
        // 1. COMPUTE k1 = f(Af)
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
        // 1.5. finally compute k1
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                k1[i*N+j][0] = agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]);
                k1[i*N+j][1] = agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]);
                
            }
        }
        
        // ------------------------------------------------------------------------------------------------------------------
        // 2. COMPUTE k2 = f(Af+k1*dt/2)
        // 2.1. assign （Af+k1*dt/2) to argument arf, compute laplacian (Af+k1*dt/2) in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Af[i*N+j][0]+k1[i*N+j][0]*dt*0.5;
                agf[i*N+j][1] = Af[i*N+j][1]+k1[i*N+j][1]*dt*0.5;
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
        // 2.5. finally compute k2
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                k2[i*N+j][0] = agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]);
                k2[i*N+j][1] = agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]);
                
            }
        }
        // ------------------------------------------------------------------------------------------------------------------
        // 3. COMPUTE k3 = f(Af+k2*dt/2)
        // 3.1. assign （Af+k2*dt/2) to argument arf, compute laplacian (Af+k2*dt/2) in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Af[i*N+j][0]+k2[i*N+j][0]*dt*0.5;
                agf[i*N+j][1] = Af[i*N+j][1]+k2[i*N+j][1]*dt*0.5;
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
        // 3.5. finally compute k3
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                k3[i*N+j][0] = agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]);
                k3[i*N+j][1] = agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]);
            }
        }
        // ------------------------------------------------------------------------------------------------------------------
        // 4. COMPUTE k4 = f(Af+k3*dt)
        // 4.1. assign （Af+k3*dt) to argument arf, compute laplacian (Af+k3*dt) in Fourier space
        for (int i=0; i<localN; ++i) {
            int globali = local0+i;
            for (int j=0; j<N; ++j) {
                agf[i*N+j][0] = Af[i*N+j][0]+k3[i*N+j][0]*dt;
                agf[i*N+j][1] = Af[i*N+j][1]+k3[i*N+j][1]*dt;
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
        // 4.5. finally compute k4
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                k4[i*N+j][0] = agf[i*N+j][0] + (2*M_PI/L)*(2*M_PI/L)*( lapf[i*N+j][0]- c1*lapf[i*N+j][1])-(nlf[i*N+j][0]+c3*nlf[i*N+j][1]);
                k4[i*N+j][1] = agf[i*N+j][1] + (2*M_PI/L)*(2*M_PI/L)*(c1*lapf[i*N+j][0] + lapf[i*N+j][1])-(nlf[i*N+j][1]-c3*nlf[i*N+j][0]);
            }
        }
        // ------------------------------------------------------------------------------------------------------------------
        // 5. update Anewf
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Anewf[i*N+j][0] = Af[i*N+j][0] + dt*1./6. *(k1[i*N+j][0] + 2.*k2[i*N+j][0] + 2.*k3[i*N+j][0] + k4[i*N+j][0]);
                Anewf[i*N+j][1] = Af[i*N+j][1] + dt*1./6. *(k1[i*N+j][1] + 2.*k2[i*N+j][1] + 2.*k3[i*N+j][1] + k4[i*N+j][1]);
            }
        }
//        for (int i=0; i<N; ++i) {
//            for (int j=0; j<N; ++j) {
//                Anewf[i*N+j][0] = Af[i*N+j][0] + dt*k1[i*N+j][0];
//                Anewf[i*N+j][1] = Af[i*N+j][1] + dt*k1[i*N+j][1];
//            }
//        }
        
        // 6. replace Af with Anewf
        for (int i=0; i<localN; ++i) {
            for (int j=0; j<N; ++j) {
                Af[i*N+j][0] = Anewf[i*N+j][0];
                Af[i*N+j][1] = Anewf[i*N+j][1];
            }
        }
        
        
    } // end of main temporal loop
    
    // get ending time
    gettimeofday(&t_end, 0);
    double t_elapsed = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)*1e-6; // compute time elapsed in seconds
    printf("The total elaspsed time is %f\n",t_elapsed);
    
    
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
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fullA = (double*)malloc(2*N*N*sizeof(double));
    }
    
    MPI_Gather(A, 2*localN*N, MPI_DOUBLE, fullA, 2*localN*N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    if (rank == 0) {
        // write sollution A to file
        FILE* fileid1 = fopen("finalsol.bin", "w");
        fwrite(fullA, sizeof(fftw_complex), N*N*2, fileid1);
        fclose(fileid1);
    }

    
    fftw_destroy_plan(pb);
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb_nl);
    fftw_destroy_plan(pf_nl);
    if (rank == 0) {
        free(fullA);
    }
    fftw_free(A);
    fftw_free(Af);
    fftw_free(Anewf);
    fftw_free(ag);
    fftw_free(agf);
    fftw_free(nl);
    fftw_free(nlf);
    fftw_free(k1);
    fftw_free(k2);
    fftw_free(k3);
    fftw_free(k4);
    MPI_Finalize();
    return 0;
}

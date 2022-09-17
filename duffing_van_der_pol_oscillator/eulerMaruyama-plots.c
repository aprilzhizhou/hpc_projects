#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char * argv[])
{
    // read in input arguments alpha, sigma, N, and M
//    double alpha = atof(argv[1]);
//    double sigma = atof(argv[2]);
//    int M = atoi(argv[3]);
//    int N = atoi(argv[4]);
    
    // calculate time - initialize time variables
    struct timeval t_start;
    struct timeval t_end;
    
    double alpha = 1;
    double sigma = 0.5;
    int N = 100000;
    int M;
    
    int Mvec[4] = {250, 500, 1000, 2000};
    
    for (int idx = 0; idx<4; ++idx)
    {
        M = Mvec[idx];
        // get starting time
        gettimeofday(&t_start, 0);
        
    // clock_t start, end;
    // double cpu_time_used;
    // start = clock();
    
    
    //printf("alpha = %f \t sigma = %f \t M = %d \t N = %d \n",alpha,sigma,M,N);
    double T = 10;
    double t;
    double dt = T/M;
    // generate two normal distributed random numbers
    long int seed = (long int)time(NULL);
    if (argc >5 )
        seed = atol(argv[5]);
    //printf("Starting seed = %ld\n",seed);
    srand48(seed);
    
    int K = M/10 + 1;
    //    printf("K = %d \n",K);
    double count[K];
    int l;
    
    // initialize count as a all zero array
    for (int i = 0; i<K;++i)
    {
        count[i] = 0;
    }
    
    
    //        // write file
    //        FILE *FILEcount1 = fopen("count1.bin","w");
    //       fwrite(count,sizeof(double),K,FILEcount1);
    //        fclose(FILEcount1);
    
    
    for (int n = 0; n<N;++n)
    {
        // define initial X and Y
        double U1 = drand48();
        double U2 = drand48();
        // define initial X and Y
        double X = sqrt(-2*log(U1))*cos(2*M_PI*U2);
        double Y = 0;
        
        if (sqrt((X+alpha)*(X+alpha)+Y*Y)<=alpha/2 || sqrt((X-alpha)*(X-alpha)+Y*Y) <=alpha/2)
        {
            count[0] = count[0]+1;
        }
        
        for (int m = 0; m<M;++m)
        {
            t = (m+1)*dt; // update time t, starting from t = 1*dt
            // create a normally distributed random number
            double V1 = drand48();
            double V2 = drand48();
            double dW = dt*(sqrt(-log(V1))*cos(2*M_PI*V2));
            X = X+Y*dt;
            Y = Y+((alpha*alpha-X*X)*X-Y)*dt+sigma*X*dW;
            
            if ( (m+1)%10 == 0 ) // count how many trials satisfies one of the inequalities
            {
                if (sqrt((X+alpha)*(X+alpha)+Y*Y)<=alpha/2 || sqrt((X-alpha)*(X-alpha)+Y*Y) <=alpha/2)
                {
                    l = (m+1)/10;
                    count[l] = count[l]+1;
                }
                
            }
        }
        
        
    }
        // get ending time
        gettimeofday(&t_end, 0);
        double t_elapsed = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)*1e-6; // compute time elapsed in seconds
        printf("%f\n",t_elapsed);
        
        // compute probabilities
    double Pvec[K];
    for (int k = 0; k<K;++k)
    {
        Pvec[k] = count[k]/N;
    }
    
    //end = clock();
    //cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;
    //printf("%f\n", cpu_time_used );
    }
    
//    // write file
//    FILE *FILEprob = fopen("Prob.out","w");
//    fwrite(Pvec,sizeof(double),K,FILEprob);
//    fclose(FILEprob);
//
    return 0;
}


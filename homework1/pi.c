#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "common.h"
#include <inttypes.h>
#include <time.h>

// #define USE_CRITICAL
#define PI 3.1415926535


void usage(int argc, char** argv);
double calcPi_Serial(int num_steps);
double calcPi_P1(int num_steps);
double calcPi_P2(int num_steps);


int main(int argc, char** argv)
{
    // get input values
    uint32_t num_steps = 100000;
    if(argc > 1) {
        num_steps = atoi(argv[1]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32"\n", num_steps);
    }
    fprintf(stdout, "The first 10 digits of Pi are %0.10f\n", PI);


    // set up timer
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();


    // calculate in serial 
    start_t = ReadTSC();
    double Pi0 = calcPi_Serial(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi serially with %"PRIu32" steps is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi0);
    
    // calculate in parallel with integration
    start_t = ReadTSC();
    double Pi1 = calcPi_P1(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" steps is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi1);


    // calculate in parallel with Monte Carlo
    start_t = ReadTSC();
    double Pi2 = calcPi_P2(num_steps);
    end_t = ReadTSC();
    printf("Time to calculate Pi in // with %"PRIu32" guesses is: %g\n",
           num_steps, ElapsedTime(end_t - start_t));
    printf("Pi is %0.10f\n", Pi2);

    
    return 0;
}


void usage(int argc, char** argv)
{
    fprintf(stdout, "usage: %s <# steps>\n", argv[0]);
}

double calcPi_Serial(int num_steps)
{
    if (num_steps <= 0) return 0.0;
    const double dx = 1.0 / (double) num_steps;
    double sum = 0.0;

    for (int i = 0; i < num_steps; ++i) {
        double x = (i + 0.5) * dx;
        sum += sqrt(1.0 - x * x);
    }
    return 4.0 * sum * dx;
}

double calcPi_P1(int num_steps)
{
    if (num_steps <= 0) return 0.0;
    const double dx = 1.0 / (double) num_steps;
    double sum = 0.0;  // shared

    #pragma omp parallel
    {
        #ifdef USE_CRITICAL
        #pragma omp for
        for (int i = 0; i < num_steps; ++i) {
            double x = (i + 0.5) * dx;
            double y = sqrt(1.0 - x*x);
            #pragma omp critical
            sum += y;
        }
        #else
        #pragma omp for
        for (int i = 0; i < num_steps; ++i) {
            double x = (i + 0.5) * dx;
            double y = sqrt(1.0 - x*x);
            #pragma omp atomic
            sum += y;
        }
        #endif
    }
    return 4.0 * sum * dx;
}

double calcPi_P2(int num_steps)
{
    if (num_steps <= 0) return 0.0;

    uint64_t inside_total = 0;

    #pragma omp parallel
    {
        unsigned int seed = (unsigned int)(time(NULL) ^ (omp_get_thread_num()*0x9e3779b9u));
        uint64_t inside = 0;

        #pragma omp for
        for (int i = 0; i < num_steps; ++i) {
            double x = (double) rand_r(&seed) / (double) RAND_MAX;
            double y = (double) rand_r(&seed) / (double) RAND_MAX;
            if (x*x + y*y <= 1.0) inside++;
        }
        #pragma omp atomic
        inside_total += inside;
    }

    return 4.0 * ((double)inside_total / (double)num_steps);
}

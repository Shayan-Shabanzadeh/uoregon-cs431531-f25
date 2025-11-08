#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "common.h"



void usage(int argc, char** argv);
void verify(int* sol, int* ans, int n);
void prefix_sum(int* src, int* prefix, int n);
void prefix_sum_p1(int* src, int* prefix, int n);
void prefix_sum_p2(int* src, int* prefix, int n);


int main(int argc, char** argv)
{
    // get inputs
    uint32_t n = 1048576;
    unsigned int seed = time(NULL);
    if(argc > 2) {
        n = atoi(argv[1]); 
        seed = atoi(argv[2]);
    } else {
        usage(argc, argv);
        printf("using %"PRIu32" elements and time as seed\n", n);

    }


    // set up data 
    int* prefix_array = (int*) AlignedMalloc(sizeof(int) * n);  
    int* input_array = (int*) AlignedMalloc(sizeof(int) * n);
    srand(seed);
    for(int i = 0; i < n; i++) {
        input_array[i] = rand() % 100;
    }


    // set up timers
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();


    // execute serial prefix sum and use it as ground truth
    start_t = ReadTSC();
    prefix_sum(input_array, prefix_array, n);
    end_t = ReadTSC();
    printf("Time to do O(N-1) prefix sum on a %"PRIu32" elements: %g (s)\n", 
           n, ElapsedTime(end_t - start_t));


    // execute parallel prefix sum which uses a NlogN algorithm
    int* input_array1 = (int*) AlignedMalloc(sizeof(int) * n);  
    int* prefix_array1 = (int*) AlignedMalloc(sizeof(int) * n);  
    memcpy(input_array1, input_array, sizeof(int) * n);
    start_t = ReadTSC();
    prefix_sum_p1(input_array1, prefix_array1, n);
    end_t = ReadTSC();
    printf("Time to do O(NlogN) //prefix sum on a %"PRIu32" elements: %g (s)\n",
           n, ElapsedTime(end_t - start_t));
    verify(prefix_array, prefix_array1, n);

    
    // execute parallel prefix sum which uses a 2(N-1) algorithm
    memcpy(input_array1, input_array, sizeof(int) * n);
    memset(prefix_array1, 0, sizeof(int) * n);
    start_t = ReadTSC();
    prefix_sum_p2(input_array1, prefix_array1, n);
    end_t = ReadTSC();
    printf("Time to do 2(N-1) //prefix sum on a %"PRIu32" elements: %g (s)\n", 
           n, ElapsedTime(end_t - start_t));
    verify(prefix_array, prefix_array1, n);


    // free memory
    AlignedFree(prefix_array);
    AlignedFree(input_array);
    AlignedFree(input_array1);
    AlignedFree(prefix_array1);


    return 0;
}

void usage(int argc, char** argv)
{
    fprintf(stderr, "usage: %s <# elements> <rand seed>\n", argv[0]);
}

void verify(int* sol, int* ans, int n)
{
    int err = 0;
    for(int i = 0; i < n; i++) {
        if(sol[i] != ans[i]) {
            err++;
        }
    }
    if(err != 0) {
        fprintf(stderr, "There was an error: %d\n", err);
    } else {
        fprintf(stdout, "Pass\n");
    }
}

void prefix_sum(int* src, int* prefix, int n)
{
    prefix[0] = src[0];
    for(int i = 1; i < n; i++) {
        prefix[i] = src[i] + prefix[i - 1];
    }
}

void prefix_sum_p1(int* src, int* prefix, int n)
{  
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        prefix[i] = src[i];
    for (int offset = 1; offset < n; offset <<= 1) {
        #pragma omp parallel for
        for (int i = (offset << 1) - 1; i < n; i += (offset << 1)) {
            prefix[i] += prefix[i - offset];
        }
    }
    prefix[n - 1] = 0;
    for (int offset = n >> 1; offset >= 1; offset >>= 1) {
        #pragma omp parallel for
        for (int i = (offset << 1) - 1; i < n; i += (offset << 1)) {
            int temp = prefix[i - offset];
            prefix[i - offset] = prefix[i];
            prefix[i] += temp;
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        prefix[i] += src[i];
}

void prefix_sum_p2(int* src, int* prefix, int n)
{
    int num_threads;
    int* partial_sums;

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int total_threads = omp_get_num_threads();

        #pragma omp single
        {
            num_threads = total_threads;
            partial_sums = malloc(sizeof(int) * num_threads);
        }

        int chunk = (n + num_threads - 1) / num_threads;
        int start = tid * chunk;
        int end = (start + chunk > n) ? n : start + chunk;

        if (start < n) {
            prefix[start] = src[start];
            for (int i = start + 1; i < end; i++)
                prefix[i] = prefix[i - 1] + src[i];

            partial_sums[tid] = prefix[end - 1];
        }

        #pragma omp barrier

        #pragma omp single
        for (int t = 1; t < num_threads; t++)
            partial_sums[t] += partial_sums[t - 1];

        #pragma omp barrier

        if (tid > 0 && start < n) {
            int offset = partial_sums[tid - 1];
            for (int i = start; i < end; i++)
                prefix[i] += offset;
        }
    }

    free(partial_sums);
}




#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define MAX_THREADS 128
#define CACHE_LINE_SIZE 64

typedef struct {
    double sum;
    char pad[CACHE_LINE_SIZE - sizeof(double)];
} padded_sum;

void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double opt_local_sum(double *x, size_t size)
{
    padded_sum local_sum[MAX_THREADS] = {0.0};
    int num_threads;

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        num_threads = omp_get_num_threads();  

        #pragma omp for
        for (size_t i = 0; i < size; i++) {
            local_sum[id].sum += x[i];
        }
    }

    double sum_val = 0.0;
    for (int i = 0; i < num_threads; i++) {
        sum_val += local_sum[i].sum;
    }
    return sum_val;
}

int main() {
    size_t size = 10000000; 
    double *input = malloc(size * sizeof(double));

    int num_trials = 10;
    double start_time, end_time, total_time;
    double times[num_trials];
    double mean_time = 0, std_dev = 0;

    // Performance measurement
    for (int i = 0; i < num_trials; i++) {
        generate_random(input, size);
        
        start_time = omp_get_wtime();
        opt_local_sum(input, size);
        end_time = omp_get_wtime();
        
        total_time = end_time - start_time;
        times[i] = total_time;
        mean_time += total_time;
    }

    mean_time /= num_trials;

    // Calculating standard deviation
    for (int i = 0; i < num_trials; i++) {
        std_dev += (times[i] - mean_time) * (times[i] - mean_time);
    }
    std_dev = sqrt(std_dev / num_trials);

    printf("Average Time: %f seconds\n", mean_time);
    printf("Standard Deviation: %f seconds\n", std_dev);

    free(input);
    return 0;
}

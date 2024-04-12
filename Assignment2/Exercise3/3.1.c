#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size)
{
  double sum_val = 0.0;

  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
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
        serial_sum(input, size);
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

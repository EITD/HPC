#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N 16 // Define the overall matrix size

void sequential_matrix_multiply(double *A, double *B, double *C)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < N; k++)
            {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

void initialize_matrix(double *mat, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            mat[i * cols + j] = i * cols + j;
        }
    }
}

void rearrange(int n, double *matrix, double *rearranged)
{
    int M = N / n; // number of block
    for (int block_i = 0; block_i < M; block_i++)
    {
        for (int block_j = 0; block_j < M; block_j++)
        { // block location
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {                                                              // element in block location
                    int src_idx = (block_i * n + i) * N + (block_j * n + j);   // element block row and col
                    int dst_idx = (block_i * M + block_j) * n * n + i * n + j; // element one-dimensional location
                    rearranged[dst_idx] = matrix[src_idx];
                }
            }
        }
    }
}

void matrix_multiply(int n, double *a, double *b, double *c)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
            {
                sum += a[i * n + k] * b[k * n + j];
            }
            c[i * n + j] += sum;
        }
    }
}

void gather(int n, double *matrix, double *rearranged)
{
    int M = N / n;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {                                                 // element location in big matrix
            int row = i / n;                              // sub-matrix row
            int col = j / n;                              // sub-matrix col
            int block = row * M + col;                    // number of sub-matrix before
            int index_in_block = (i % n) * n + (j % n);   // element offset in sub-matrix
            int index = block * (n * n) + index_in_block; // element location in one-dimensional
            rearranged[i * N + j] = matrix[index];
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, size, dims[2], periods[2], reorder, coords[2], provided;
    int up, down, left, right;
    MPI_Comm cart_comm;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    dims[0] = dims[1] = sqrt(size);
    int n = N / dims[0];
    if (N != n * dims[0])
    {
        if (rank == 0)
        {
            printf("The square root of the number of processes required by this program should be divisible by N.\n");
            printf("The current number of processes is: %d, The square root of the number of processes is: %d, N is: %d\n", size, dims[0], N);
        }
        MPI_Finalize();
        return 0;
    }
    periods[0] = periods[1] = 1;
    reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    double *A = calloc(n * n, sizeof(double));
    double *B = calloc(n * n, sizeof(double));
    double *C = calloc(n * n, sizeof(double));
    double *tempA = calloc(n * n, sizeof(double));
    double *full_A = NULL;
    double *full_B = NULL;
    double *full_C = NULL;
    double *rearrange_full_A = NULL;
    double *rearrange_full_B = NULL;
    double *rearrange_full_C = NULL;
    double *C_seq = NULL;

    if (rank == 0)
    {
        full_A = malloc(N * N * sizeof(double));
        full_B = malloc(N * N * sizeof(double));
        full_C = calloc(N * N, sizeof(double));
        initialize_matrix(full_A, N, N);
        initialize_matrix(full_B, N, N);

        // normal matrix multiply
        C_seq = malloc(N * N * sizeof(double));
        sequential_matrix_multiply(full_A, full_B, C_seq);

        rearrange_full_A = malloc(N * N * sizeof(double));
        rearrange_full_B = malloc(N * N * sizeof(double));
        rearrange(n, full_A, rearrange_full_A);
        rearrange(n, full_B, rearrange_full_B);
    }

    // int grid_rows = N / n;
    // int row_index = rank / grid_rows;
    // int col_index = rank % grid_rows;
    // int starts[2] = {row_index * n, col_index * n};

    // MPI_Datatype blocktype;
    // int sizes[2] = {N, N};  // Size of full matrix
    // int subsizes[2] = {n, n};  // Size of submatrix
    // int starts[2] = {0, 0};  // Starting coordinates for each submatrix
    // MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &blocktype);
    // MPI_Type_create_resized(blocktype, 0, n*sizeof(double), &blocktype);
    // MPI_Type_commit(&blocktype);

    // int *sendcounts = malloc(size * sizeof(int));
    // int *displs = malloc(size * sizeof(int));

    // for (int i = 0; i < size; i++) {
    //     sendcounts[i] = 1;
    //     int row_block = i / (N / n);
    //     int col_block = i % (N / n);
    //     displs[i] = row_block * n * N + col_block * n;
    // }

    // int sendcounts[4] = {1, 1, 1, 1};
    // int displs[4] = {0, 2, 8, 10};

    // MPI_Scatterv(full_A, sendcounts, displs, blocktype, A, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Scatterv(full_B, sendcounts, displs, blocktype, B, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Scatter(rearrange_full_A, n * n, MPI_DOUBLE, A, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(rearrange_full_B, n * n, MPI_DOUBLE, B, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Comm row_comm;
    MPI_Comm_split(cart_comm, coords[0], coords[1], &row_comm);

    // Perform the algorithm
    for (int k = 0; k < dims[0]; k++)
    {
        int bcast_root = (coords[0] + k) % dims[0]; // Root shifts rightward each step
        if (coords[1] == bcast_root)
        {
            memcpy(tempA, A, n * n * sizeof(double));
        }
        MPI_Bcast(tempA, n * n, MPI_DOUBLE, bcast_root, row_comm);
        matrix_multiply(n, tempA, B, C);

        // Roll the B matrix upwards
        MPI_Cart_shift(cart_comm, 0, -1, &down, &up);
        MPI_Sendrecv_replace(B, n * n, MPI_DOUBLE, up, 0, down, 0, cart_comm, MPI_STATUS_IGNORE);
    }

    MPI_Gather(C, n * n, MPI_DOUBLE, full_C, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        rearrange_full_C = malloc(N * N * sizeof(double));
        gather(n, full_C, rearrange_full_C);

        printf("Fox result:\n");
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                printf("%6.2f ", rearrange_full_C[i * N + j]);
            }
            printf("\n");
        }

        printf("-------------------------------------------------------------------------\n");

        printf("Sequential result:\n");
        int count = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                printf("%6.2f ", C_seq[i * N + j]);
                if (rearrange_full_C[i * N + j] != C_seq[i * N + j])
                {
                    count++;
                }
            }
            printf("\n");
        }

        printf("-------------------------------------------------------------------------\n");
        printf("Different count: %d\n", count);

        free(C_seq);
        free(full_A);
        free(full_B);
        free(full_C);
        free(rearrange_full_A);
        free(rearrange_full_B);
        free(rearrange_full_C);
    }

    free(A);
    free(B);
    free(C);
    free(tempA);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
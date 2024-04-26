#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    int rank, size, i, provided;

    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain


    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells
    int nxn_loc = nxc/size + 2; // number of nodes is number cells + 1; we add also 2 ghost cells
    // printf("%d %d %d\n", nxn_loc, nxc, size);
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);

    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);

    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
    double left_ghost, right_ghost;
    MPI_Request send_request, recv_request;
    MPI_Status status;

    if (size == 1) {
      right_ghost = f[1];
      left_ghost = f[nxn_loc - 2];
    } else {
      if (rank == 0) {
        MPI_Isend(&f[nxn_loc - 2], 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &send_request);
        MPI_Irecv(&left_ghost, 1, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Wait(&recv_request, &status);
        MPI_Wait(&send_request, &status);

        MPI_Isend(&f[1], 1, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &send_request);
        MPI_Irecv(&right_ghost, 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Wait(&recv_request, &status);
        MPI_Wait(&send_request, &status);
      } else {
        MPI_Irecv(&left_ghost, 1, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Isend(&f[nxn_loc - 2], 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &send_request);
        MPI_Wait(&recv_request, &status);
        MPI_Wait(&send_request, &status);

        MPI_Irecv(&right_ghost, 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Isend(&f[1], 1, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &send_request);
        MPI_Wait(&recv_request, &status);
        MPI_Wait(&send_request, &status);
      }
    }

    f[0] = left_ghost;
    f[nxn_loc - 1] = right_ghost;
    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    // Print f values
    if (rank==0 || rank==1 || rank == (size-1)) { // print only rank 0 for convenience
      printf("My rank %d of %d\n", rank, size );
      printf("Here are f including ghost cells (my value, dfdx, cos)\n");
      for (i=0; i<nxn_loc; i++) {
        if (i == 1 || i == (nxn_loc-1)) {
          printf("--------------------\n");
        }
        printf("%f, %f, %f\n", f[i], dfdx[i], cos(L_loc*rank + (i-1) * dx));
      }
      printf("\n");
    }

    free(f);
    free(dfdx);

    MPI_Finalize();
}


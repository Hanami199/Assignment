#include <stdio.h>
#include <mpi.h>


int main( int argc, char **argv) 
{   
    int Ntasks, Myrank;
    int provided;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &Myrank);

    int N = 0;

    if (Myrank == 0) {
        N = 3;
        MPI_Send( &N, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if (Myrank == 1){
        printf("%d\n", N);
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("%d\n", N);
    }

    MPI_Finalize();
    return 0;
}
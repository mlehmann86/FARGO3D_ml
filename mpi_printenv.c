#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // For gethostname

// Use extern char **environ; (common, works on most Linux)
extern char **environ;

int main(int argc, char *argv[]) {
    int rank, size;
    char hostname[256];
    char **env; // Pointer for iterating

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    gethostname(hostname, sizeof(hostname));

    // Add a barrier to try and sync output slightly
    MPI_Barrier(MPI_COMM_WORLD);

    // Have only rank 0 print a header
    if (rank == 0) {
        printf("--- Environment variables visible to MPI processes ---\n");
        printf("--- Format: [Rank / Size @ Hostname] Variable=Value ---\n\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Each rank iterates through and prints its environment
    env = environ;
    while (*env) {
        // Print format: [Rank/Size @ Hostname] Variable=Value
        printf("[%d/%d @ %s] %s\n", rank, size, hostname, *env);
        env++;
    }

    MPI_Barrier(MPI_COMM_WORLD);
     if (rank == 0) {
        printf("\n--- End of Environment Dump ---\n");
    }

    MPI_Finalize();
    return 0;
}

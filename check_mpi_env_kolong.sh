#!/bin/bash

# Script to compile and run a small MPI program to dump environment variables
# on Kolong to find the local rank variable for Fargo3D.

echo "--- Starting MPI Environment Check for Kolong ---"
START_DIR=$(pwd)
TMP_DIR=$(mktemp -d -t mpi_env_check_XXXXXX)

if [ ! -d "$TMP_DIR" ]; then
    echo "ERROR: Could not create temporary directory!"
    exit 1
fi

cd "$TMP_DIR" || exit 1
echo ">>> Working in temporary directory: ${TMP_DIR} <<<"

# 1. Define C source code using a here document
echo ">>> Creating mpi_printenv.c <<<"
cat << 'EOF' > mpi_printenv.c
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
EOF

if [ ! -f mpi_printenv.c ]; then
  echo "ERROR: Failed to create mpi_printenv.c in ${TMP_DIR}"
  cd "$START_DIR" || exit 1
  rm -rf "$TMP_DIR"
  exit 1
fi
echo "mpi_printenv.c created successfully."
echo ""

# 2. Load necessary modules
echo ">>> Loading modules (cuda/11.7, compiler/2022.1.0, mpi/latest) <<<"
module purge || echo "Warning: module purge failed, continuing..."
module add cuda/11.7
if [ $? -ne 0 ]; then echo "ERROR: Failed to load cuda/11.7"; cd "$START_DIR" || exit 1; rm -rf "$TMP_DIR"; exit 1; fi
module add compiler/2022.1.0
if [ $? -ne 0 ]; then echo "ERROR: Failed to load compiler/2022.1.0"; cd "$START_DIR" || exit 1; rm -rf "$TMP_DIR"; exit 1; fi
module add mpi/latest
if [ $? -ne 0 ]; then echo "ERROR: Failed to load mpi/latest"; cd "$START_DIR" || exit 1; rm -rf "$TMP_DIR"; exit 1; fi
echo "Modules loaded."
module list # Show loaded modules
echo ">>> Setting I_MPI_CC=icc <<<"
export I_MPI_CC=icc
echo ""

# 3. Compile the C code
echo ">>> Compiling mpi_printenv.c with mpicc <<<"
mpicc mpi_printenv.c -o mpi_printenv
if [ $? -ne 0 ]; then
  echo "ERROR: Compilation failed."
  cd "$START_DIR" || exit 1
  rm -rf "$TMP_DIR" # Clean up
  exit 1
fi
if [ ! -f mpi_printenv ]; then
    echo "ERROR: Compiled executable mpi_printenv not found!"
    cd "$START_DIR" || exit 1
    rm -rf "$TMP_DIR" # Clean up
    exit 1
fi
echo "Compilation successful: mpi_printenv executable created."
echo ""

# 4. Run the compiled program with mpirun
OUTPUT_FILE="${START_DIR}/kolong_env_output.txt" # Save to original directory
echo ">>> Running with 'mpirun -np 2 ./mpi_printenv' <<<"
echo "Output will be saved to ${OUTPUT_FILE}"
mpirun -np 2 ./mpi_printenv > "${OUTPUT_FILE}"
if [ $? -ne 0 ]; then
  echo "ERROR: mpirun execution failed."
  # Keep output file for inspection
  cd "$START_DIR" || exit 1
  rm -rf "$TMP_DIR" # Clean up source and exe
  exit 1
fi
echo "Execution finished."
echo ""

# 5. Output guidance and potential candidates
echo "--- Analysis ---"
if [ ! -s "${OUTPUT_FILE}" ]; then
    echo "ERROR: Output file ${OUTPUT_FILE} is empty or missing!"
else
    echo "Full output saved to ${OUTPUT_FILE}"
    echo "Please examine ${OUTPUT_FILE} carefully."
    echo "Look for variables containing 'RANK', 'LOCAL', 'PMI', 'MPI', 'HYDRA', 'SLURM'."
    echo "The LOCAL RANK variable should have value 0 for one rank and 1 for the other (if on the same node)."
    echo ""
    echo ">>> Potential candidate variables found by grep (check values in file!): <<<"
    # Grep for relevant patterns, ignore common irrelevant ones, show unique lines prefixed by rank
    # Sort primarily by variable name after '=' sign, then by rank
    grep -E 'RANK|LOCAL|PMI|MPI|HYDRA|SLURM' "${OUTPUT_FILE}" \
      | grep -vE 'MPI_SUFFIX|MODULE|MANPATH|INFO|LIB|INC|DIR|PATH|PS1|CPU_BIND|PIN|DOMAIN' \
      | sed -e 's/\[\([0-9]*\)\/\([0-9]*\) @ \(.*\)\] \(.*\)/\1 \4/' \
      | sort -k2 -k1n \
      | uniq
    echo "----------------------------------------------------------------------------"
fi
echo ""

# 6. Cleanup
echo ">>> Cleaning up temporary directory (${TMP_DIR}) <<<"
cd "$START_DIR" || exit 1
rm -rf "$TMP_DIR"
echo ""

echo "--- MPI Environment Check Done ---"

#!/bin/bash

# Path to the executable
EXEC="../bin/parallel_solver_mpi"

# Number of dimensions
NUM_DIM=2

# Number of processes
NUM_PROCESSES=4

# Grid size
M_DIM=16
N_DIM=$M_DIM
K_DIM=$M_DIM

# Number of processors in each dimension
M_PROCESSORS=2
N_PROCESSORS=2
K_PROCESSORS=1

# Change to the parent directory
cd ..

# Clean previous build
make clean

# Compile the program
make all

# Change back to the original directory
cd tests

if [ $? -eq 0 ]; then
  # Compilation successful, run the program
  mpirun -np $NUM_PROCESSES $EXEC $NUM_DIM $M_DIM $N_DIM $M_PROCESSORS $N_PROCESSORS
else
  # Compilation failed
  echo "Compilation failed."
fi


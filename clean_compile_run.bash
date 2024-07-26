#!/bin/bash

# Path to the executable
EXEC="exec/parallel_solver_mpi.out"

# Number of processes
NUM_PROCESSES=8

# Clean output
rm -f output/*
rm -f python/output/*

# Clean previous build
make clean

# Compile the program either in debug or release mode
make release

if [ $? -eq 0 ]; then
  # Compilation successful, run the program
  mpirun --report-bindings --bind-to socket -np $NUM_PROCESSES $EXEC
else
  # Compilation failed
  echo "Compilation failed."
fi

# jupyter nbconvert --execute --inplace python/read_system_solution.ipynb


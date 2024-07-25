# Program name
PROGNAME = parallel_solver_mpi

# Compiler
FC = mpifort
COMPILERSTANDARD = #-std=gnu

# Source and binary directories
SRC_DIR = src
BIN_DIR = bin
INC_DIR = include
EXEC_DIR = exec

### Compilation flags ###

# Debug flags
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fbounds-check -O0

# Release flags
RELEASEFLAGS = -O3 -march=native -flto -funroll-loops

# Verbose optimization flags
VERBOSEFLAGS = -fopt-info-vec -fopt-info-loop -fopt-info-inline

# OpenMP flags
OPENMPFLAGS = -fopenmp

# ThreadSanitizer flag
SANITIZERFLAGS = -fsanitize=thread

# Other flags
OTHERFLAGS = -c
EFLAGS = -fimplicit-none -Wall -Wextra #-Werror
IFLAGS = -J$(INC_DIR)/
MFLAGS = -m64 -fPIC -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

# Linker flags for linking
LDFLAGS = -lm -lopenblas -lscalapack $(OPENMPFLAGS) #$(SANITIZERFLAGS)

# Executable name
EXEC = $(EXEC_DIR)/$(PROGNAME).out

# Files to compile
FILES = utility_functions_module.f90 constants_module.f90 mpi_wrapper_module.f90 functions_module.f90 comm_module.f90 block_module.f90 initialization_module.f90 FD_module.f90 multigrid_module.f90 solver_module.f90 FD_test_module.f90 block_test_module.f90 test_Poisson_module.f90 nonlinear_test_module.f90 Lid_driven_cavity_benchmark.f90 TravellingWave_2D.f90 scalapack_module.f90 main.f90

# Source files with directory prefix
SRCS = $(addprefix $(SRC_DIR)/,$(FILES))
# Object files
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BIN_DIR)/%.o,$(SRCS))

# Targets
.PHONY: default release clean

default:
	@echo "To compile the program type make debug|release"

# Combined compiler flags for compilation
release: CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(RELEASEFLAGS) $(EFLAGS) $(IFLAGS) $(OPENMPFLAGS) $(MFLAGS) $(VERBOSEFLAGS) #$(SANITIZERFLAGS)
release: $(EXEC)

debug: CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(DEBUGFLAGS) $(EFLAGS) $(IFLAGS) $(MFLAGS) $(VERBOSEFLAGS) #$(SANITIZERFLAGS)
debug: $(EXEC)

$(EXEC): $(OBJS) | $(EXEC_DIR)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR) $(INC_DIR)
	$(FC) $(CFLAGS) $< -o $@

$(BIN_DIR) $(INC_DIR) $(EXEC_DIR):
	mkdir -p $@	

clean:
	rm -f $(EXEC)
	rm -f $(BIN_DIR)/*.o 
	rm -f $(INC_DIR)/*.mod 
	rm -f $(SRC_DIR)/*.mod

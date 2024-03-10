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
RELEASEFLAGS = -O3 -flto

# OpenACC flags
OPENACCFLAGS = -fopenacc

# Other flags
OTHERFLAGS = -c
EFLAGS = -fimplicit-none -Wall -Wextra #-Werror
IFLAGS = -J$(INC_DIR)/
MFLAGS = -m64 -fPIC -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

# Linker flags for linking
LDFLAGS = -llapack -lblas

# Executable name
EXEC = $(EXEC_DIR)/$(PROGNAME).out

# Files to compile
FILES = utility_functions_module.f90 mpi_wrapper_module.f90 constants_module.f90 functions_module.f90 comm_module.f90 block_module.f90 initialization_module.f90 finite_difference_module.f90 rank_module.f90  solver_module.f90 main.f90

# Source files with directory prefix
SRCS = $(addprefix $(SRC_DIR)/,$(FILES))
# Object files
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BIN_DIR)/%.o,$(SRCS))

# Targets
.PHONY: default release clean

default:
	@echo "To compile the program type make debug|release"

# Combined compiler flags for compilation
release: CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(RELEASEFLAGS) $(EFLAGS) $(IFLAGS) $(OPENACCFLAGS) $(MFLAGS)
release: $(EXEC)

debug: CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(DEBUGFLAGS) $(EFLAGS) $(IFLAGS) $(OPENACCFLAGS) $(MFLAGS)
debug: $(EXEC)

$(EXEC): $(OBJS) | $(EXEC_DIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR) $(INC_DIR)
	$(FC) $(CFLAGS) $< -o $@

$(BIN_DIR) $(INC_DIR) $(EXEC_DIR):
	mkdir -p $@	

clean:
	rm -f $(EXEC)
	rm -f $(BIN_DIR)/*.o 
	rm -f $(INC_DIR)/*.mod 
	rm -f $(SRC_DIR)/*.mod

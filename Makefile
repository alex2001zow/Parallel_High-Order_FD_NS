# Compiler
FC = mpifort
COMPILERSTANDARD = #-std=gnu

# Source and binary directories
SRC_DIR = src
BIN_DIR = bin
INC_DIR = include
EXEC_DIR = exec

# Compilation flags
OTHERFLAGS = -c
OFLAGS = -g -O0
EFLAGS = -fimplicit-none -Wall -Wextra #-Werror
IFLAGS = -J$(INC_DIR)/
MFLAGS = -m64 -fopenmp -fPIC -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

# Combined compiler flags for compilation
CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(OFLAGS) $(EFLAGS) $(IFLAGS) $(MFLAGS)

# Linker flags for linking
LDFLAGS = 

# Executable name
EXEC = $(EXEC_DIR)/parallel_solver_mpi.out

# Files to compile
FILES = utility_functions_module.f90 mpi_wrapper_module.f90 constants_module.f90 neighbor_types_module.f90 comm_module.f90 block_module.f90 initialization_module.f90 rank_module.f90 solver_module.f90 main.f90

# Source files with directory prefix
SRCS = $(addprefix $(SRC_DIR)/,$(FILES))
# Object files
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BIN_DIR)/%.o,$(SRCS))

# Targets
.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS) | $(EXEC_DIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR) $(INC_DIR)
	$(FC) $(CFLAGS) $< -o $@

$(BIN_DIR) $(INC_DIR) $(EXEC_DIR):
	mkdir -p $@

clean:
	rm -f $(EXEC) $(BIN_DIR)/*.o $(INC_DIR)/*.mod $(SRC_DIR)/*.mod

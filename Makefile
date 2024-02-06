# Compiler
FC = mpifort
COMPILERSTANDARD = #-std=gnu

# Source and binary directories
SRC_DIR = src
BIN_DIR = bin
INC_DIR = include

# Compilation flags
OTHERFLAGS = -c
OFLAGS = -O3
EFLAGS = -fimplicit-none -Wall -Wextra #-Werror
IFLAGS = -J$(INC_DIR)/
MFLAGS = -m64 -fopenmp -fPIC -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

# Combined compiler flags for compilation
CFLAGS = $(COMPILERSTANDARD) $(OTHERFLAGS) $(OFLAGS) $(EFLAGS) $(IFLAGS) $(MFLAGS)

# Linker flags for linking
LDFLAGS = 

# Executable name
EXEC = $(BIN_DIR)/parallel_solver_mpi.out

# Files to compile
FILES = mpi_wrapper_module.f90 constants_module.f90 neighbor_types_module.f90 utility_functions_module.f90 initilization_module.f90 rank_parameters_module.f90 solver_module.f90 print_module.f90 main.f90

# Source files with directory prefix
SRCS = $(addprefix $(SRC_DIR)/,$(FILES))
# Object files
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BIN_DIR)/%.o,$(SRCS))

# Targets
.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR) $(INC_DIR)
	$(FC) $(CFLAGS) $< -o $@

$(BIN_DIR) $(INC_DIR):
	mkdir -p $@

clean:
	rm -f $(EXEC) $(BIN_DIR)/*.o $(INC_DIR)/*.mod $(SRC_DIR)/*.mod
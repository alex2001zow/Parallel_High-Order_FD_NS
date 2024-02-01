# Compiler
FC = mpifort

# Compilation flags
OTHERFLAGS = -c

# Optimization flags
OFLAGS = -O3

# Error flags
EFLAGS = -Wall -Wextra

# Include flags
IFLAGS = -Jinclude/

# Module flags
MFLAGS = -m64 -fopenmp -fPIC -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8

# Combined compiler flags for compilation
CFLAGS = $(OTHERFLAGS) $(OFLAGS) $(EFLAGS) $(IFLAGS) $(MFLAGS)

# Linker flags for linking
LDFLAGS = 

# Executable name
EXEC = bin/parallel_solver_mpi

# Source files
SRCS = src/mpi_wrapper_module.f90 src/constants_module.f90 src/utility_functions_module.f90 src/rank_parameters_module.f90 src/main.f90

# Object files
OBJS = $(patsubst src/%.f90,bin/%.o,$(SRCS))

# Targets
all: $(EXEC)

bin/:
	mkdir -p bin/

include/:
	mkdir -p include/

$(EXEC): bin/ $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

bin/%.o: src/%.f90 include/
	$(FC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) bin/*.o include/*.mod src/*.mod
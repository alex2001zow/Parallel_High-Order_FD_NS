# Compiler
FC = mpifort

# Compiler flags for compilation
CFLAGS = -c -O3 -Wall -Wextra -Jinclude/

# Linker flags for linking
LDFLAGS = 

# Executable name
EXEC = bin/parallel_solver_mpi

# Source files
SRCS = src/constants_module.f90 src/utility_functions_module.f90 src/rank_parameters_module.f90 src/main.f90

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
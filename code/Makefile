# Compiler
FC = mpifort

# Compiler flags for compilation
CFLAGS = -c -O3 -Wall -Wextra

# Linker flags for linking
LDFLAGS = 

# Executable name
EXEC = parallel_solver_mpi

# Source files
SRCS = rank_parameters_module.f90 utility_functions_module.f90 main.f90

# Object files
OBJS = $(patsubst %.f90,%.o,$(SRCS))

# Targets
all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)

%.o: %.f90
	$(FC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC) $(OBJS)

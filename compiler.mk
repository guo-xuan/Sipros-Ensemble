# GNU Compiler for OpenMP installations
CC := g++
OPTS := -std=c++11 -fopenmp -O3 -g

# Intel Compiler for OpenMP installations
# CC := icpc
# OPTS := -std=c++11 -fopenmp -O3 -g

# MPI C++ Compiler for OpenMPI installations (default)
# MCC := mpic++
# MOPTS := -std=c++11 -fopenmp -O3 -g

# Intel MPI C++ Compiler installations
# MCC := icpc
# MOPTS := -std=c++11 -openmp -O3 -g -lmpi -lmpi++

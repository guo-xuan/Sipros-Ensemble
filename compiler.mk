# GNU Compiler for OpenMP installations
CC := g++
GCC = gcc
OPTS := -std=c++11 -fopenmp -O3 -g

# Intel Compiler for OpenMP installations
# CC := icpc
# OPTS := -std=c++11 -fopenmp -O3 -g

# MPI C++ Compiler for OpenMPI installations (default)
MCC := /home/xgo/local/bin/mpic++
MGCC := /home/xgo/local/bin/mpicc
MOPTS := -std=c++11 -fopenmp -O3 -g

# Intel MPI C++ Compiler installations
# MCC := icpc
# MGCC := icc
# MOPTS := -std=c++11 -qopenmp -O3 -g -lmpi -lmpi++

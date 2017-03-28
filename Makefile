#############################################
### MAKE file for Sipros Ensemble 
#############################################

OpenMP_DIR = ReleaseOpenMP/
MPI_DIR = ReleaseMPI/

openmp:
	$(MAKE) -C $(OpenMP_DIR)
	mkdir -p bin
	cp $(OpenMP_DIR)/Sipros_OpenMP bin/
	
mpi:
	$(MAKE) -C $(MPI_DIR)
	mkdir -p bin
	cp $(MPI_DIR)/Sipros_MPI bin/

all:
	$(MAKE) -C $(OpenMP_DIR)
	$(MAKE) -C $(MPI_DIR)
	mkdir -p bin
	cp $(OpenMP_DIR)/Sipros_OpenMP bin/
	cp $(MPI_DIR)/Sipros_MPI bin/
	
clean:
	$(MAKE) -C $(OpenMP_DIR) clean
	$(MAKE) -C $(MPI_DIR) clean
	-$(RM) -rf bin/Sipros_OpenMP bin/Sipros_MPI

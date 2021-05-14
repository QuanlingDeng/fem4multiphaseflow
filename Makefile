TRGTS = LINALG MULTIFRONTAL UMFPACK MESH FEM GENERAL TIME MSROBIN
MPI_TRGTS = MPI_LINALG MPI_MESH MPI_FEM 

all: $(TRGTS)

mpi_all: $(MPI_TRGTS)

LINALG:
	cd linalg ; make

MULTIFRONTAL:
	cd linalg/multifrontal ; make

UMFPACK:
	cd linalg/multifrontal/UMFPACK2.2 ; make

MESH:
	cd mesh ; make

FEM:
	cd fem ; make

MSROBIN:
	cd msrobin ; make

GENERAL:
	cd general ; make

TIME:
	cd timedependent ; make

MPI_LINALG:
	cd linalg ; make mpi_all

MPI_MESH:
	cd mesh ; make mpi_all

MPI_FEM:
	cd fem ; make mpi_all

veryclean:
	cd linalg ; make veryclean
	cd mesh ; make veryclean
	cd fem ; make veryclean 
	cd msrobin ; make veryclean 	
	cd general ; make veryclean
	cd linalg/multifrontal ; make veryclean
	cd linalg/multifrontal/UMFPACK2.2 ; make veryclean
	cd timedependent ; make veryclean

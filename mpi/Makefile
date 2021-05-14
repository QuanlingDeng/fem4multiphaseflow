include ../comp.mk
mpi_objects = mpi_lib.o

%.o: %.cpp
	$(MPICC) -c -o $@ $< 

all:	$(mpi_objects)
	$(AR) -rcs lib_mpi.a $(mpi_objects)	

veryclean:
	rm -f *.exe *.o *.a *~
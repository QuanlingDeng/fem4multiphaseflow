include ../comp.mk
mesh_objects = point.o segment.o triangle.o quadrilateral.o \
               tetrahedral.o vertex.o mesh.o dualmesh.o \
               controlvolume.o \

mpi_mesh_objects = point.o segment.o triangle.o quadrilateral.o \
               tetrahedral.o vertex.o mesh_p.o

%_p.o: %_p.cpp
	$(MPICC) -c -o $@ $< 

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<  

all: $(mesh_objects)
	$(AR) -rcs libmesh.a $(mesh_objects)

mpi_all: $(mpi_mesh_objects)
	$(AR) -rcs mpi_lib_mesh.a $(mpi_mesh_objects)	

veryclean:
	rm -f *.exe *.o *.a *~

include ../comp.mk
fem_objects = Data.o fem.o postprocessing.o basisfunctions.o \
	      linearfem1d.o quadraticfem1d.o cubicfem1d.o \
              linearfem2d.o quadraticfem2d.o cubicfem2d.o bilinearfem2d.o linearfemelast2d.o bilinearfemelast2d.o \
	      linearfvem2d.o quadraticfvem2d.o cubicfvem2d.o linearfemgeomech2d.o \
	      linearfem3d.o \
	      #stokes2dMINI.o stokes2dP1P0.o stokes2dCR.o quadraticbsplinefem2d.o 

mpi_fem_objects = Data.o fem_p.o linearfem1d_p.o quadraticfem1d_p.o cubicfem1d_p.o \
                  linearfem2d_p.o linearfem3d_p.o bilinearfem2d_p.o quadraticfem2d.o \

%_p.o: %_p.cpp
	$(MPICC) -c -o $@ $< 

%.o: %.cpp
	$(CC) $(CFLAGS) $(GSL) -c  -o $@ $<  

all: $(fem_objects)
	$(AR) -rcs libfem.a $(fem_objects)

mpi_all:$(mpi_fem_objects)
	$(AR) -rcs mpi_lib_fem.a $(mpi_fem_objects)		

veryclean:
	rm -f *.exe *.o *.a *~

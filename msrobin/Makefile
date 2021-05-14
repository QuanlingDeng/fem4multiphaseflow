include ../comp.mk
msrobin_objects = msbdry.o msd.o directMsRM.o \

mpi_msrobin_objects = msd.o msbdry.o \

%_p.o: %_p.cpp
	$(MPICC) -c -o $@ $< 

%.o: %.cpp
	$(CC) $(CFLAGS) -c  -o $@ $<  

all: $(msrobin_objects)
	$(AR) -rcs libmsrobin.a $(msrobin_objects)

mpi_all:$(mpi_msrobin_objects)
	$(AR) -rcs mpi_lib_msrobin.a $(mpi_msrobin_objects)		

veryclean:
	rm -f *.exe *.o *.a *~

include ../comp.mk

linalg_objects = vector.o \
             matrix.o \
             tridiagmat.o \
             rectmat.o \
             cgsolver.o \
             pcgsolver.o \
             bicgstab.o \
             gmres.o \
             sparsemat.o \
             sparsesmoothers.o \
             inversemat.o \
             symdensemat.o \
             bandmat.o

%.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $< 

all:	$(linalg_objects)
	$(AR) -rcs lib_linalg.a $(linalg_objects)

veryclean:
	rm -f *.exe *.o *.a *~


include ../comp.mk
timedependent_objects = hyperbolicproblem.o  firstorderupwind.o secondorderupwind.o parabolicproblem.o fdata.o twophaseflow.o \

%.o: %.cpp
	$(CC) $(CFLAGS) -c  -o $@ $<  

all: $(timedependent_objects)
	$(AR) -rcs libtimedependent.a $(timedependent_objects)		

veryclean:
	rm -f *.exe *.o *.a *~

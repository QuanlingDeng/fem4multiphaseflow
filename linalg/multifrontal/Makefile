include ../../comp.mk

multiobj = multifrontal.o

%.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $< 

all:	$(multiobj)
	$(AR) -rcs libmultifrontal.a $(multiobj)

veryclean:
	rm -f *.exe *.o *.a *~




include ../comp.mk
mesh_objects = container.o patch.o stack.o geometry.o 
%.o: %.cpp
	$(CC) $(CFLAGS) $(GSL) -c -o $@ $<  

all: $(mesh_objects)
	$(AR) -rcs  libgeometry.a  $(mesh_objects) 	

veryclean:
	rm -f *.exe *.o *.a *~


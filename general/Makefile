# C++ compiler (g++,CC), options
include ../Compiler.mk

# general code library source files
LIBCPPS = array.cpp \
	tic_toc.cpp \
	table.cpp \
	stable3d.cpp \
	sort_pairs.cpp \
	osockstream.cpp \
	mergesort.cpp \
	sets.cpp \
	evalfunc.cpp function.cpp trigaussquad.cpp

LIBOBJS = $(LIBCPPS:.cpp=.o)

# Library will be linked with -l$(LNAME)
LNAME   = general

# Library rules
include ../Library.mk

# C++ compiler (g++,CC)
CPP = g++

MYHOMEDIR	= /home/qdeng

# GNU project C++ 
CPPOPTSg++   = -Wno-deprecated -Wall
OPTOPTSg++   = -O3
DEBUGOPTSg++ = -g -DDEBUG -fno-inline
# for some Sun machines
LINKOPTSg++  = -z muldefs 
# for isc sgi machines
LINKOPTSg++  = -Wl,-woff,131
# lechery, hilbert, fourier
LINKOPTSg++  = 

# MIPSpro C++ 
CPPOPTSCC    = -w -ptused -lm
OPTOPTSCC    = -O3 -lm
DEBUGOPTSCC  = -g -DDEBUG -lm

# Portland Group's C++ compiler
CPPOPTSpgCC   = -tused
OPTOPTSpgCC   = +K3 -O -fast
DEBUGOPTSpgCC = -g -DDEBUG
LINKOPTSpgCC  = -w --implicit_include --one_instantiation_per_object -tused


CPPOPTS   = $(CPPOPTS$(CPP)) 
OPTOPTS   = $(OPTOPTS$(CPP))
LINKOPTS  = $(LINKOPTS$(CPP))

# Rule for compiling source files
.SUFFIXES: .cpp 
.cpp.o:
	$(CPP) $(CPPOPTS) -c ${INCLUDE} $< 

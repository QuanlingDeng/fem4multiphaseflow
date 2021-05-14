# Library names
LIBNAME    = lib$(LNAME).a

# ar and ranlib format
AROPTS        = cruv
AR_NO_RANLIB  = @echo  "Building $@..."; \
                ar $(AROPTS) $(LIBNAME) $(LIBOBJS)
AR_AND_RANLIB = @echo  "Building $@..."; \
                ar $(AROPTS) $(LIBNAME) $(LIBOBJS); \
	        ranlib $(LIBNAME)

# List different architectures
ARiris4d     = $(AR_NO_RANLIB)
ARcray       = $(AR_NO_RANLIB)

ARsun4       = $(AR_AND_RANLIB)
ARi586-linux = $(AR_AND_RANLIB)
ARi386-linux = $(AR_AND_RANLIB)
ARi386       = $(AR_AND_RANLIB)
ARalpha      = $(AR_AND_RANLIB)
ARrs6000     = $(AR_AND_RANLIB)
ARi686       = $(AR_AND_RANLIB)
ARx86_64     = $(AR_AND_RANLIB)
ARx86_64-linux = $(AR_AND_RANLIB)
ARiris4d     = $(AR_NO_RANLIB)
ARcray       = $(AR_NO_RANLIB)
ARintel-pc = $(AR_AND_RANLIB)


# Use $(HOSTTYPE) to choose the architecture
AR = $(ARx86_64)


###########################  TARGETS  #############################

# default target
all       : $(LIBNAME)


# Build the library
$(LIBNAME): $(LIBOBJS) 
	$(AR)

deps:
	$(CPP) $(CPPOPTS) -MM $(LIBCPPS) > deps.mk

# Clean unnecessary files
clean:
	rm -f core *.o

# Clean also the libraries, backups, etc
veryclean:
	rm -f core *.o *.a *~
	rm -rf ii_files ti_files SunWS_cache

ALL = adcp_Linux-x86_64

OS = $(shell uname -s)
CFLAGS = -std=c99 -O2 # -D_GNU_SOURCE #-fgnu89-inline
CPP=g++
CPPFLAGS = -Wall -O2
OPENMPFLAGS = -fopenmp
LDFLAGS = -lm
LDFLAGS_DEBUG = -lm

ifeq ($(OS), Linux)
	CFLAGS = -std=c99 -Wall -O2 #-Wno-unused-result 
	CFLAGS_DEBUG = -std=c99 -Wall -O0 -g -DDEBUG #-Wno-unused-result -g
	ifneq ($(shell which mpicc),)
		MPICC = mpicc
		MPILDFLAGS = $(LDFLAGS)
		ALL := $(ALL) peptmpi
	endif
endif

ifeq ($(OS), Darwin)
	CFLAGS = -std=c99 -Wall -O2
	CFLAGS_DEBUG = -std=c99 -Wall -O0 -g -arch i386 -DDEBUG
	ifneq ($(shell which mpicc),)
		MPICC = mpicc
		MPILDFLAGS = $(LDFLAGS)
		ALL := $(ALL) peptmpi
	endif
endif

ifeq ($(OS), SunOS)
	CFLAGS = -xO2
	ifneq ($(shell which mpcc),)
		MPICC = mpcc
		MPILDFLAGS = -lmpi $(LDFLAGS)
		ALL := $(ALL) peptmpi
	endif
endif

all : $(ALL)

#serial peptide program (MC, nested sampling)
adcp_Linux-x86_64 : nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c canonicalAA.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@ -g

clean :
	$(RM) $(ALL) $(TOOLS)

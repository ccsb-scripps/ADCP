ALL = rama coma oops stats befa lipa peptide merg dssp2cm peptmpi cdlearn
TOOLS = tools/extract_distances \
      tools/1pga_helix_orientation \
      tools/R0003_collect_SS \
      tools/num_helix_residues \
      tools/count_charged_residues \
      tools/icm2icm \
      tools/initialize_test \
      tools/turn_analyser \
      tools/hbond_pattern \
      tools/rmsd/rmsd \
      tools/rmsd/rmsd_debug \
      tools/rmsd/rmsd_dihedral \
      tools/energy_landscape_charts/list_large_basins \
      tools/checkpoint_split/checkpoint_split \
      tools/thermodynamics/thermo \
      tools/energy_landscape_charts/beta \
      tools/energy_landscape_charts/rmsdbasin \
      tools/energy_landscape_charts/rmsdbasin_dihedral

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

tools : $(TOOLS)

rama : aadict.c ramachandran.c vector.c peptide.c rotation.c params.c error.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

coma : aadict.c cm.c params.h error.h 
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

oops : oops.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

stats : statistics.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

befa : aadict.c bfactor.c peptide.c rotation.c vector.c params.c error.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

lipa : aadict.c pauling.c peptide.c rotation.c vector.c params.c error.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

#internal cdlearn program
cdlearn : aadict.c energy.c cdlearn.c metropolis.c flex.c peptide.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DLJ_HBONDED_HARD -DLJ_NEIGHBOUR_HARD $(OPENMPFLAGS) $^ $(LDFLAGS) -o $@

cdlearn_debug : aadict.c energy.c cdlearn.c metropolis.c flex.c peptide.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS_DEBUG) -DLJ_HBONDED_HARD -DLJ_NEIGHBOUR_HARD $^ $(LDFLAGS_DEBUG) -o $@

#serial peptide program (MC, nested sampling)
peptide : nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c canonicalAA.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@ -g

peptide_debug : nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS_DEBUG) $^ $(LDFLAGS_DEBUG) -o $@

#parallel peptmpi program (parallel tempering, nested sampling)
peptmpi : nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c random16.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(MPICC) $(CFLAGS) -DPARALLEL $^ $(MPILDFLAGS) -o $@

peptmpi_debug : nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c random16.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(MPICC) $(CFLAGS_DEBUG) -DPARALLEL $^ $(LDFLAGS_DEBUG) -o $@

merg : mergie.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

dssp2cm : dssp2cm.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

#TOOLS

tools/extract_distances : aadict.c energy.c tools/tools_main.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DEXTRACT_DISTANCES $^ $(LDFLAGS) -o $@

tools/1pga_helix_orientation : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DHELIX_ORIENTATION $^ $(LDFLAGS) -o $@

tools/R0003_collect_SS : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DR0003_SSBOND $^ $(LDFLAGS) -o $@

tools/num_helix_residues : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c probe.c flex.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DNUM_HELIX_RESIDUES $^ $(LDFLAGS) -o $@

tools/count_charged_residues : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DCOUNT_CHARGED_RESIDUES $^ $(LDFLAGS) -o $@

tools/icm2icm : tools/icm2icm.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

tools/initialize_test : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DINITIALIZE_TEST $^ $(LDFLAGS) -o $@

tools/turn_analyser : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c flex.c probe.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DTURN_ANALYSER $^ $(LDFLAGS) -o $@

tools/hbond_pattern : tools/tools_main.c aadict.c energy.c metropolis.c peptide.c probe.c flex.c rotation.c vector.c params.c error.c checkpoint_io.c vdw.c
	$(CC) $(CFLAGS) -DHBOND_PATTERN $^ $(LDFLAGS) -o $@

tools/rmsd/rmsd : aadict.c peptide.c tools/rmsd/rmsd.c tools/rmsd/rmsd_main.c rotation.c vector.c flex.c params.c metropolis.c energy.c error.c vdw.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

tools/rmsd/rmsd_debug : aadict.c peptide.c tools/rmsd/rmsd.c tools/rmsd/rmsd_main.c rotation.c vector.c flex.c params.c metropolis.c energy.c error.c vdw.c
	$(CC) $(CFLAGS_DEBUG) $^ $(LDFLAGS_DEBUG) -o $@

tools/rmsd/rmsd_dihedral : aadict.c peptide.c tools/rmsd/rmsd.c tools/rmsd/rmsd_main.c flex.c rotation.c vector.c params.c metropolis.c energy.c error.c vdw.c
	$(CC) $(CFLAGS) -DRMS_DIHEDRAL $^ $(LDFLAGS) -o $@

tools/energy_landscape_charts/list_large_basins : tools/energy_landscape_charts/list_large_basins.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

tools/checkpoint_split/checkpoint_split : tools/checkpoint_split/checkpoint_split.c tools/checkpoint_split/rmsd.c params.c vector.c rotation.c aadict.c peptide.c metropolis.c energy.c error.c vdw.c
	$(CC) $(CFLAGS)  $^ $(LDFLAGS) -o $@

#CPP TOOLS

tools/thermodynamics/thermo : tools/thermodynamics/thermo.cpp
	$(CPP) $(CPPFLAGS) $^ $(LDFLAGS) -o $@

tools/energy_landscape_charts/beta : tools/energy_landscape_charts/beta.cpp
	$(CPP) $(CPPFLAGS) $^ $(LDFLAGS) -o $@

tools/energy_landscape_charts/rmsdbasin : tools/energy_landscape_charts/rmsdbasin.cpp
	$(CPP) $(CPPFLAGS) $(OPENMPFLAGS) $^ $(LDFLAGS) -o $@

tools/energy_landscape_charts/rmsdbasin_dihedral : tools/energy_landscape_charts/rmsdbasin.cpp
	$(CPP) $(CPPFLAGS) $(OPENMPFLAGS) -DRMS_DIHEDRAL $^ $(LDFLAGS) -o $@

wrap :
	tar -zcvf crankite.tar.gz README ChangeLog TODO COPYING Makefile *.[ch] viewer.ps *.sh *.awk engh-huber-2001.txt peptide-contacts.txt *.py tests tools

clean :
	$(RM) $(ALL) $(TOOLS)

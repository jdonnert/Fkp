## P-Fkp Makefile ## 
SHELL = /bin/bash

# Target Computer
ifndef SYSTYPE
	SYSTYPE	:= $(shell hostname)
endif

#Std systype
CC		 	= mpicc
OPTIMIZE 	= -Wall -g -O3
MPI_LIBS 	=	-lmpi
MPI_INCL 	= 
GSL_INCL 	=
GSL_LIBS 	= 
PGPLOT_INCL = 
PGPLOT_LIBS =

ifeq ($(SYSTYPE),geeuw.strw.leidenuniv.nl)
CC           =  mpicc
OPTIMIZE     = -O3 -m64 -Wall -g -march=native -mtune=native -flto -fwhole-program
MPI_LIBS     = -lmpich -lrt 
MPI_INCL     = 
GSL_INCL     = -I/home/jdonnert/Libs/include
GSL_LIBS     = -L/home/jdonnert/Libs/lib
PGPLOT_INCL  =
PGPLOT_LIBS  =
endif

ifeq ($(SYSTYPE),MSI)
CC           =  mpicc
OPTIMIZE     = -O2 -m64 -Wall -g -xHost
MPI_LIBS     = -lmpich -lrt -L /home/jonestw/donne219/Libs/$(MYMACHINE)/lib
MPI_INCL     = -I /home/jonestw/donne219/Libs/$(MYMACHINE)/include
GSL_INCL     =
GSL_LIBS     =
PGPLOT_INCL  =
PGPLOT_LIBS  =
endif

ifeq ($(SYSTYPE),para33.strw.leidenuniv.nl)
	CC           =  mpicc
OPTIMIZE     = -O3  -m64 -g -march=native -flto -fwhole-program
	MPI_LIBS     = -lmpich -lrt
	MPI_INCL     =
	GSL_INCL     =
	GSL_LIBS     =
	PGPLOT_INCL  =
	PGPLOT_LIBS  =
endif


ifeq ($(SYSTYPE),mach64.ira.inaf.it)
CC           =  mpicc
OPTIMIZE     = -O3  -m64 -g -march=bdver1 -mtune=native -mprefer-avx128 -mieee-fp -minline-all-stringops -fprefetch-loop-arrays --param prefetch-latency=300 -funroll-all-loops -flto -fwhole-program
MPI_LIBS     =  -lmpich
MPI_INCL     =
GSL_INCL     = -I/home/donnert/Libs/include
GSL_LIBS     = -L/home/donnert/Libs/lib
PGPLOT_INCL  =
PGPLOT_LIBS  =
endif



ifeq ($(SYSTYPE),"DARWIN")
CC           =  mpicc
OPTIMIZE     = -O3 -std=c99 -mtune=native -march=corei7
MPI_LIBS     = -lmpich -L/Users/julius/Devel/lib
MPI_INCL     = -I/Users/julius/Devel/include
GSL_INCL     = 
GSL_LIBS     = 
PGPLOT_INCL  =
PGPLOT_LIBS  =
endif

ifeq ($(SYSTYPE),"MPA")
CC           =  mpicc
OPTIMIZE     = -O3 -Wall -g
MPI_LIBS     = -mpi 
MPI_INCL     = 
GSL_INCL     = -I/afs/mpa/home/jdonnert/Libs/@sys/include
GSL_LIBS     = -L/afs/mpa/home/jdonnert/Libs/@sys/lib
PGPLOT_INCL  = -I/usr/local/include
PGPLOT_LIBS  = -L/usr/common/pdsoft/appl/pgplot/ -L/usr/X11R6/lib
endif

ifeq ($(SYSTYPE),"RZG-ODIN") 
CC       =  mpicc
OPTIMIZE = -O3  -Wall -g
MPI_LIBS = 
MPI_INCL = 
GSL_INCL = -I/u/jdonnert/Libs/include
GSL_LIBS = -L/u/jdonnert/Libs/lib
endif


# Target
EXEC 	= P-Fkp		

# Sources 
SRCDIR	= src/

OBJFILES= io/io.o \
		  		io/txt.o \
				io/gadget.o \
				io/binary.o \
		  modules/zero.o	\
				modules/secondaries_brunetti.o \
				modules/simple_secondaries.o \
				modules/hard_sphere.o \
				modules/brunetti_07.o \
				modules/donnert_13.o \
				modules/expansion.o \
				modules/shock_primaries.o \
		 common.o main.o setup.o solver.o cosmo.o timestep.o sort.o  \
		 sort_particles.o timing.o compress.o unit.o print_settings.o
OBJS	= $(addprefix $(SRCDIR),$(OBJFILES))

INCLFILES	= common.h cosmo.h compress.h unit.h constants.h config.h sort.h \
			  timing.h io/io.h \
			  		io/gadget.h \
			  modules/modules.h	\
			  ../Makefile
INCL	= $(addprefix $(SRCDIR),$(INCLFILES))

LIBS = -lm -lgsl -lgslcblas  $(MPI_LIBS) $(GSL_LIBS) 
CFLAGS = -std=c99  -fopenmp $(OPTIMIZE)


$(EXEC)	: $(OBJS)
	$(CC) $(CFLAGS)  $(OBJS)  $(LIBS) -o $(EXEC)
	ctags -R src/*.[ch]

$(OBJS)	: $(INCL)

$(SRCDIR)config.h : Config 
	sed '/^#/d; /^$$/d; s/^/#define /g' Config > $(SRCDIR)config.h

.ONESHELL:
$(SRCDIR)print_settings.c : Config
	echo '#include "common.h"' >  $(SRCDIR)print_settings.c
	echo 'void print_compile_time_settings(){' >> $(SRCDIR)print_settings.c
	echo 'rprintf("Compiled with : \n"' >> $(SRCDIR)print_settings.c
	sed '/^#/d; /^$$/d; s/^/"   /g; s/$$/ \\n"/g;' Config >>  $(SRCDIR)print_settings.c
	echo ');}' >> $(SRCDIR)print_settings.c

clean	: 
	rm -f  $(OBJS) $(EXEC)



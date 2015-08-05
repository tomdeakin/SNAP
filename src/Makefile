#
# Build SNAP
#

TARGET	= snap

ifndef COMPILER
	MESSAGE=select a compiler to compiler in OpenMP, e.g. make COMPILER=INTEL
endif


OMP_INTEL	= -openmp
OMP_CRAY	=
OMP_GNU		= -fopenmp

OMP=$(OMP_$(COMPILER))

FFLAGS_		= -O3
FFLAGS_CRAY	= -O3
FFLAGS_INTEL	= -O3
FFLAGS_GNU	= -O3

CFLAGS_		= -O3
CFLAGS_CRAY	= -O3
CFLAGS_INTEL	= -O3 -std=c99
CFLAGS_GNU	= -O3 -std=c99

MPI_COMPILER	= mpif90
C_MPI_COMPILER	= mpicc


FFLAGS=$(FFLAGS_$(COMPILER)) $(OMP) $(OPTIONS)
CFLAGS=$(CFLAGS_$(COMPILER)) $(OMP) $(OPTIONS)

.SUFFIXES:	.f90 .o

OBJS = global.o snap_main.o utils.o version.o plib.o geom.o sn.o \
       data.o control.o input.o setup.o dealloc.o translv.o solvar.o \
       outer.o expxs.o inner.o sweep.o octsweep.o dim1_sweep.o dim3_sweep.o \
       output.o time.o mms.o omp_sweep.o omp_sweep_c.o

SRCS = global.f90 snap_main.f90 utils.f90 version.f90 plib.f90 geom.f90 \
       sn.f90 data.f90 control.f90 input.f90 setup.f90 dealloc.f90 \
       translv.f90 solvar.f90 outer.f90 expxs.f90 inner.f90 sweep.f90 \
       octsweep.f90 dim1_sweep.f90 dim3_sweep.f90 output.f90 time.f90 mms.f90 \
       omp_sweep.f90

$(TARGET) :	$(OBJS)
		$(MPI_COMPILER) $(FFLAGS) -o $@ $(OBJS)
		@echo $(MESSAGE);

include make.deps

#
# Fortran rule
#
%.o:	%.f90 Makefile
	$(MPI_COMPILER) $(FFLAGS) -c $<

%.o:	%.c Makefile
	$(C_MPI_COMPILER) $(CFLAGS) -c $<

#
# Cleanup
#
clean:
	rm -f *.o *.mod *.bc

#
# Count lines of code
#
count:
	rm -f Lines
	for file in $(SRCS); do ./LineCount $$file Lines; done
	gawk -f ./LineReport < Lines >> Lines
	cat Lines

#
# Link compiled files only. No recompile.
#
link:
	$(FORTRAN) $(FFLAGS) $(FFLAG2) -o $(TARGET) $(OBJS)



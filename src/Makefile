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
		echo $(MESSAGE);

snap_main.o :	global.o utils.o version.o input.o setup.o output.o time.o \
		dealloc.o plib.o control.o omp_sweep.o
utils.o :	global.o dealloc.o plib.o
version.o :	global.o
plib.o :	global.o time.o
geom.o :	global.o
sn.o :		global.o
data.o :	global.o
control.o :	global.o
input.o :	global.o plib.o geom.o sn.o data.o control.o utils.o time.o
setup.o :	global.o geom.o sn.o data.o control.o mms.o utils.o time.o \
		plib.o
dealloc.o :	global.o sn.o data.o mms.o geom.o solvar.o control.o
translv.o :	global.o plib.o geom.o sn.o data.o control.o utils.o \
		solvar.o expxs.o outer.o time.o
solvar.o :	global.o plib.o geom.o sn.o data.o control.o
outer.o :	global.o geom.o sn.o data.o solvar.o control.o inner.o \
		time.o plib.o expxs.o
expxs.o :	global.o geom.o sn.o data.o
inner.o :	global.o geom.o sn.o data.o control.o solvar.o sweep.o \
		time.o plib.o
sweep.o :	global.o geom.o sn.o data.o control.o octsweep.o solvar.o
octsweep.o :	global.o geom.o sn.o data.o control.o solvar.o \
		dim1_sweep.o dim3_sweep.o
dim1_sweep.o :	global.o geom.o sn.o data.o control.o
dim3_sweep.o :	global.o plib.o geom.o sn.o data.o control.o
output.o :	global.o geom.o sn.o data.o control.o utils.o solvar.o \
		mms.o time.o plib.o
time.o :	global.o
mms.o:	global.o geom.o sn.o data.o control.o plib.o

omp_sweep.o:	global.o geom.o data.o sn.o omp_sweep.f90 omp_sweep_c.o


#
# Fortran rule
#
%.o:	%.f90
	$(MPI_COMPILER) $(FFLAGS) -c $<

%.o:	%.c
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



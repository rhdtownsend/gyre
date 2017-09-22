# File     : Make.inc
# Purpose  : standard makefile rules, targets

# The following flags can be edited (set to either yes or no)

# Enable debugging (with a performance penalty)
DEBUG?=no

# Enable OpenMP parallelization
OMP?=yes

# Enable MPI parallelization
MPI?=no

# Use double precision arithmetic (best left to 'yes')
DOUBLE_PRECISION?=yes

# Use IEEE fortran features
IEEE?=no

# Enable FPE checks
FPE?=yes

# Enable experimental features
EXPERIMENTAL?=no

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# Preprocesser flags

OPTS=DEBUG OMP MPI HDF5 DOUBLE_PRECISION IEEE

FPX3FLAGS:=${FPX3FLAGS} -I.:${SUBDIRS} $(foreach opt,${OPTS},$(if $(filter yes,${${opt}}),-D${opt}))

# Compiler settings

export F9XC

ifeq (${MPI},yes)
  F9XC=mpif90
  LDFLAGS:=${LDFLAGS} -lpciaccess -lnuma
else
  F9XC=gfortran
endif

MODPATH=${MESASDK_ROOT}/include

FFLAGS:=${FFLAGS} -finit-real=snan -fbacktrace -fmax-errors=25

ifeq (${FPE},yes)
  FFLAGS:=${FFLAGS} -ffpe-trap=invalid,overflow,zero
endif

ifeq (${DEBUG},yes)
  FFLAGS:=${FFLAGS} -fcheck=all -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized -finline-limit=0 -ggdb
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=DEBUG -DGFORTRAN_PR_69185 -DGFORTRAN_PR_49636
else
  FFLAGS:=${FFLAGS} -O2 -march=native
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=INFO -DGFORTRAN_PR_49636
endif

F9XFLAGS:=${F9XFLAGS} $(addprefix -I,${MODPATH}) ${FFLAGS} -std=f2008
F77FLAGS:=${F77FLAGS} ${FFLAGS} -ffixed-form

ifeq (${OMP},yes)
  F9XFLAGS:=${F9XFLAGS} -fopenmp
  F77FLAGS:=${F77FLAGS} -fopenmp
endif

LDFLAGS:=${LDFLAGS} `mesasdk_lapack95_link`

# Rules

vpath %.mod ${MODPATH}
vpath %.fpp ${SUBDIRS}
vpath %.inc ${SUBDIRS}
vpath %.f ${SUBDIRS}

.PRECIOUS : %.f90

%.o : %.mod

%.f90 : %.fpp
	@fpx3 ${FPX3FLAGS} < $< > $@

%.o : %.f
	@echo FC $<
	@$(F9XC) ${F77FLAGS} -c $<

%.o : %.f90
	@echo FC $<
	@${F9XC} ${F9XFLAGS} ${MODINCS} -c $<

%.mod : %.o
	@true

%.a :
	@echo AR $@
	@rm -f $@
	@${AR} crs $@ $^

% : %.o
	@echo LD $@
	@${F9XC} ${F9XFLAGS} -o $@ $^ ${LDFLAGS}

all : ${TARGETS}

install : ${TARGETS}
	@for t in ${TARGETS}; do \
            echo CP $${t} ${BINDIR}; \
            cp $${t} ${BINDIR}; \
        done

clean :
	rm -f ${TARGETS} $(addsuffix .f90,${TARGETS}) *.o *.mod .depend *.f90

almostclean :
	rm -f ${TARGETS} $(addsuffix .f90,${TARGETS}) *.o *.mod *.f90

# Dependencies

ifndef DEPEND_FILES
  DEPEND_FILES:=${TARGETS}
endif

ifneq ($(MAKECMDGOALS),clean)
-include .depend
endif

.depend :
	@echo Determining dependencies
	@fpx3_deps ${FPX3FLAGS} ${DEPEND_FILES} > .depend

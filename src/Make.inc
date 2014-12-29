# File     : Make.inc
# Purpose  : standard makefile rules, targets

# The following flags can be edited (set to either yes or no)

# Enable debugging (with a performance penalty)
DEBUG=no

# Enable OpenMP parallelization
OMP=yes

# Enable MPI parallelization
MPI=no

# Use double precision arithmetic (best left to 'yes')
DOUBLE_PRECISION=yes

# Use IEEE fortran features
IEEE=no

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
else
  F9XC=gfortran
endif

MODPATH=${MESASDK_ROOT}/include

FFLAGS:=${FFLAGS} -finit-real=snan -ffpe-trap=invalid,overflow,zero -fbacktrace 

ifeq (${DEBUG},yes)
  FFLAGS:=${FFLAGS} -fcheck=all -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized -finline-limit=0 -g
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=DEBUG
else
  FFLAGS:=${FFLAGS} -O2 -march=native -mno-avx -mno-bmi2
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=INFO
endif

F9XFLAGS:=${F9XFLAGS} $(addprefix -I,${MODPATH}) ${FFLAGS} -std=f2008
F77FLAGS:=${F77FLAGS} ${FFLAGS} -ffixed-form

ifeq (${OMP},yes)
  F9XFLAGS:=${F9XFLAGS} -fopenmp
  F77FLAGS:=${F77FLAGS} -fopenmp
endif

GFORTRAN_MAJOR:=$(shell gfortran -dumpversion | awk '{split($$0,verarray,".")} END{print verarray[1]}')
GFORTRAN_MINOR:=$(shell gfortran -dumpversion | awk '{split($$0,verarray,".")} END{print verarray[2]}')

ifneq (${GFORTRAN_MAJOR},4)
  VERSION_WARN:='y'
else
  ifneq (${GFORTRAN_MINOR},10)
    VERSION_WARN:='y'
  endif
endif

ifeq (${VERSION_WARN},y)
  $(warning You need gfortran 4.10.0 to compile this code)
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

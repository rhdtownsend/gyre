# $Id: Make.inc 446 2011-10-06 21:59:22Z townsend $

# Preprocesser flags

OPTS=DEBUG MPI OMP DOUBLE_PRECISION TIMER

FPX3FLAGS:=${FPX3FLAGS} -I.:${SUBDIRS} $(foreach opt,${OPTS},$(if ${${opt}},-D${opt}))

# Compiler settings

export F9XC

ifdef MPI
  F9XC=mpif90
else
  F9XC=gfortran
endif

MODPATH=${MESASDK_ROOT}/include

FFLAGS:=-finit-real=snan -ffpe-trap=invalid,overflow,zero -fbacktrace

ifdef DEBUG
  FFLAGS:=${FFLAGS} -fcheck=all -Wall -finline-limit=0 -g
else
  FFLAGS:=${FFLAGS} -O2 -march=native
endif

F9XFLAGS:=${F9XFLAGS} $(addprefix -I,${MODPATH}) ${FFLAGS} -std=f2008
F77FLAGS:=${F77FLAGS} ${FFLAGS} -ffixed-form

ifdef OMP
F9XFLAGS:=${F9XFLAGS} -fopenmp
F77FLAGS:=${F77FLAGS} -fopenmp
endif

LDFLAGS:=${LDFLAGS}
FPX3FLAGS:=${FPX3FLAGS} -DGFORTRAN_PR_56218

# Rules

vpath %.mod ${MODPATH}
vpath %.fpp ${SUBDIRS}
vpath %.f90 ${SUBDIRS}
vpath %.inc ${SUBDIRS}
vpath %.f ${SUBDIRS}

.PRECIOUS : %.f90

%.o : %.mod

%.f90 : %.fpp
	fpx3 ${FPX3FLAGS} < $< > $@

%.o %.mod : %.f90
	${F9XC} ${F9XFLAGS} ${MODINCS} -c $<
	@if [ -e $(basename $@).mod ]; then touch $(basename $@).mod; fi

%.o : %.f
	$(F9XC) ${F77FLAGS} -c $<

% : %.o
	${F9XC} ${F9XFLAGS} -o $@ $^ ${LDFLAGS}

all : ${TARGETS}

clean :
	rm -f ${TARGETS} $(addsuffix .f90,${TARGETS}) *.o *.mod \
          .depend *.f90
	for d in $(subst :, ,${SUBDIRS}); do if [ -e $$d/Makefile ]; then make -C $$d clean; fi; done
	

# Dependencies

ifneq ($(MAKECMDGOALS),clean)
-include .depend
endif

.depend :
	fpx3_deps ${FPX3FLAGS} ${TARGETS} > .depend

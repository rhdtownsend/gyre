# File     : Make.inc
# Purpose  : standard makefile rules, targets

# The following flags can be edited (set to either yes or no)

# Enable debugging (with a performance penalty)
DEBUG?=no

# Enable OpenMP parallelization
OMP?=yes

# Enable MPI parallelization
MPI?=no

# Enable correctly-rounded math
CRMATH?=yes

# Enable portable math (for bit-for-bit reproducibility; setting to
# no may give a performance boost)
PORTABLE?=yes

# Use IEEE fortran features
IEEE?=yes

# Enable FPE checks
FPE?=yes

# Enable experimental features
EXPERIMENTAL?=no

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# Check GYRE_DIR is set

ifneq ($(MAKECMDGOALS),clean)
  ifeq (${GYRE_DIR},) 
    $(error Please set the GYRE_DIR environment variable to point to the top-level GYRE directory)
  endif
endif

# Check SDK version

MINIMUM_VERSION = 23.7.2

ifneq ($(MAKECMDGOALS),clean)

  VERSION_CHECK = $(shell ./check_sdk_version $(MINIMUM_VERSION))

  ifneq ($(VERSION_CHECK),passed)
    ifeq ($(VERSION_CHECK),SDK too old)
      $(error Your MESA SDK installation is too old; please update to version $(MINIMUM_VERSION) or newer)
    else ifeq ($(VERSION_CHECK),SDK not installed)
      $(error MESA SDK is either not installed or not properly initialized)
    else
      $(error Could not determine MESA SDK availability)
    endif
  endif

endif

# Make flags

SHELL := /bin/bash

# Preprocessor flags

OPTS = DEBUG OMP FPE IEEE

FYPPFLAGS := ${FYPPFLAGS} \
             $(foreach subdir,$(subst :, ,${SUBDIRS}),-I${subdir}) \
	     $(foreach opt,${OPTS},$(if $(filter yes,${${opt}}),-D${opt}=1,-D${opt}=0))

# Compiler settings

# Fortran

export F9XC

F9XC = gfortran

FFLAGS := ${FFLAGS} -finit-real=snan -fbacktrace -fmax-errors=25 -fPIC

ifeq (${FPE},yes)
  FFLAGS := ${FFLAGS} -ffpe-trap=invalid,overflow,zero
endif

ifeq (${DEBUG},yes)
  FFLAGS := ${FFLAGS} -Og -ggdb -fcheck=all -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized -finline-limit=0
  FYPPFLAGS := ${FYPPFLAGS} -DLOG_LEVEL=\'DEBUG\' -DGFORTRAN_PR_69185 -DGFORTRAN_PR_49636
else
  FFLAGS := ${FFLAGS} -O2 -march=native
  FYPPFLAGS := ${FYPPFLAGS} -DLOG_LEVEL=\'INFO\' -DGFORTRAN_PR_49636
endif

F9XFLAGS := ${F9XFLAGS} $(addprefix -I,${MODPATH}) ${FFLAGS} -std=f2018
F77FLAGS := ${F77FLAGS} ${FFLAGS} -ffixed-form

ifeq (${OMP},yes)
  F9XFLAGS := ${F9XFLAGS} -fopenmp
  F77FLAGS := ${F77FLAGS} -fopenmp
endif

ifeq ($(PORTABLE),yes)
  F9XFLAGS := ${F9XFLAGS} -ffp-contract=off
  F77FLAGS := ${F9XFLAGS} -ffp-contract=off
else
  F9XFLAGS := ${F9XFLAGS} -march=native
  F77FLAGS := ${F77FLAGS} -march=native
endif

ifeq ($(CRMATH),yes)
  export MESASDK_MATH_SLOT=crmath
  LDFLAGS := ${LDFLAGS} `crmath_link` `crlibm_link`
else
  export MESASDK_MATH_SLOT=default
endif

# Linking

LDFLAGS := ${LDFLAGS} `hdf5_link` `lapack95_link`

ifeq ($(shell uname -s),Darwin)
  LIB_SUFFIX = .dylib
  LIB_LDFLAGS = -shared -install_name ${GYRE_DIR}/lib/$@
  EXEC_LDFLAGS = 
else			   
  LIB_SUFFIX = .so
  LIB_LDFLAGS = -Wl,--no-warn-execstack -shared -Wl,-rpath -Wl,${GYRE_DIR}/lib
  EXEC_LDFLAGS = -Wl,--no-warn-execstack -Wl,-rpath -Wl,${GYRE_DIR}/lib
endif

# Targets-related variables

LIB_TARGETS_ = $(addsuffix ${LIB_SUFFIX},${LIB_TARGETS})

DEP_TARGETS = $(patsubst %,.%.d,${LIB_TARGETS} ${EXEC_TARGETS})

MOD_TARGETS = $(patsubst lib%,%_m.mod,${LIB_TARGETS})
INC_TARGETS = $(patsubst lib%,%.h,$(filter libc%,${LIB_TARGETS}))

# Paths

vpath %.fypp ${SUBDIRS}
vpath %.inc ${SUBDIRS}
vpath %.pyx ${SUBDIRS}
vpath %.c ${SUBDIRS}
vpath %.h ../src/include

# Explicit rules

.PHONY: all install lib-install exec-install headers-install clean almostclean

all : ${LIB_TARGETS_} ${EXEC_TARGETS}

install : lib-install exec-install headers-install

lib-install : ${LIB_TARGETS_}
	@for t in $^; do \
            echo CP $${t} ${LIBDIR}; \
            cp $${t} ${LIBDIR}; \
        done

exec-install : ${EXEC_TARGETS}
	@for t in $^; do \
            echo CP $${t} ${BINDIR}; \
            cp $${t} ${BINDIR}; \
        done

headers-install : ${MOD_TARGETS} ${INC_TARGETS}
	@for t in $^; do \
            echo CP $${t} ${INCDIR}; \
            cp $${t} ${INCDIR}; \
        done

clean : almostclean
	rm -f .*.d .*.dpp

almostclean :
	rm -f ${LIB_TARGETS_} ${EXEC_TARGETS} *.o *.mod *.smod *.f90 *.c

${LIB_TARGETS_} :
	@echo LD_LIB $@
	@${F9XC} ${F9XFLAGS} -o $@ $^ ${LDFLAGS} ${LIB_LDFLAGS}

${EXEC_TARGETS} :
	@echo LD $@
	@${F9XC} ${F9XFLAGS} -o $@ $^ ${LDFLAGS} ${EXEC_LDFLAGS}

${DEP_TARGETS} : 
	@echo DEPEND $(patsubst .%.d,%,$@) "->" $@
	@makedepf90 $(addprefix -d ,${EXTRA_DEPS}) -m %m.mod -Wconfused $^ > $@
	@echo "include \$$(patsubst %.f90,.%.dpp,$^)" >> $@

# Implicit rules

%.o : %.mod
%.o : %.c
% : %.mod
% : %.o

%.f90 .%.dpp : %.fypp
	@./fypp_deps ${FYPPFLAGS} $< $*.f90 > .$*.dpp

%_m.o %_m.mod : %_m.f90
	@echo FC $<
	@${F9XC} ${F9XFLAGS} -c $<

%.o : %.f90
	@echo FC $<
	@${F9XC} ${F9XFLAGS} -c $<

%.d : EXTRA_DEPS = $(addsuffix .mod,${$*_DEPS})

.PRECIOUS : %.f90 %.c

# Macros

strip_d = $(patsubst %.d,%,$(1))

# Dependency handling

ifneq ($(MAKECMDGOALS),clean)
-include ${DEP_TARGETS}
endif

.SECONDEXPANSION:

${LIB_TARGETS_} : %${LIB_SUFFIX} : $${$$*_OBJS} $$(addsuffix ${LIB_SUFFIX},$${$$*_LIBS})

${EXEC_TARGETS} : $$($$@_OBJS) $$(addsuffix ${LIB_SUFFIX},$${$$@_LIBS})

${DEP_TARGETS} : .%.d : $$(subst .o,.f90,$${$$*_OBJS})

# File     : Make.inc
# Purpose  : standard makefile rules, targets

############
# Variables
############

SRC_DIRS = $(addprefix ../src/, ad bvp cad common context core diff	\
   ext extern frontend/gyre grid interp math math/crmath matrix mode	\
   model nad output par poly rad rot sad search tar tide tnad)

DEPEND_DIR = .depend

ALL_TARGETS = $(EXEC_TARGETS) $(LIB_TARGETS)

ALL_SRCS = $(foreach target,$(ALL_TARGETS),$($(target)_SRCS))

ALL_DPPS = $(addprefix $(DEPEND_DIR)/,$(ALL_SRCS:%.fypp=%.dpp))
ALL_DEPS = $(addprefix $(DEPEND_DIR)/,$(ALL_TARGETS:%=%.dep))

##########
# Targets
##########

all: $(ALL_TARGETS)

$(DEPEND_DIR) :
	mkdir -p $(DEPEND_DIR)

clean :
	rm -f $(ALL_TARGETS) $(DEPEND_DIR)/* *.f90 *.mod *.smod *.o 

.PHONY: clean

.SUFFIXES:

########
# Paths
########

vpath %.fypp $(SRC_DIRS)

#################
# Preproccessing
#################

# Flags & files

FPP = ./fypp_deps

FPPFLAGS := $(FPPFLAGS) \
            $(addprefix -I,$(SRC_DIRS)) \
            $(foreach opt,DEBUG OMP FPE IEEE,$(if $(filter yes,$($(opt))),-D$(opt)=1,-D$(opt)=0))

ifeq ($(DEBUG),yes)
  FPPFLAGS := $(FPPFLAGS) -DLOG_LEVEL=\'DEBUG\' -DGFORTRAN_PR_69185 -DGFORTRAN_PR_49636
else
  FPPFLAGS := $(FPPFLAGS) -DLOG_LEVEL=\'INFO\' -DGFORTRAN_PR_49636
endif

# Rules

%.f90 : %.fypp $(DEPEND_DIR)/%.dpp | $(DEPEND_DIR)
	@$(FPP) $(FPPFLAGS) -d $(DEPEND_DIR)/$*.dpp $< $@

$(ALL_DPPS) : 

# Extra dependencies

ifneq ($(MAKECMDGOALS),clean)
  -include $(ALL_DPPS)
endif

############
# Compiling
############

# Flags & files

FC = gfortran
MAKEDEP = awk -f makedepf08.awk

FFLAGS := $(FFLAGS) $(addprefix -I,$(MODPATH)) -std=f2018 \
          -finit-real=snan -fbacktrace -fmax-errors=25 -fPIC

ifeq ($(FPE),yes)
  FFLAGS := $(FFLAGS) -ffpe-trap=invalid,overflow,zero
endif

ifeq ($(DEBUG),yes)
  FFLAGS := $(FFLAGS) -Og -ggdb -fcheck=all -Wall \
            -Wno-unused-dummy-argument -Wno-maybe-uninitialized -finline-limit=0
else
  FFLAGS := $(FFLAGS) -O2
endif

ifeq ($(OMP),yes)
  FFLAGS := $(FFLAGS) -fopenmp
endif

ifeq ($(PORTABLE),yes)
  FFLAGS := $(FFLAGS) -ffp-contract=off
else
  FFLAGS := $(FFLAGS) -march=native
endif

# Rules

%.o : $(DEPEND_DIR)/%.anc
	@echo FC $*
	@$(FC)$ $(FFLAGS) -c $*.f90

$(DEPEND_DIR)/%.anc : %.f90 | $(DEPEND_DIR)
	@$(FC) $(FFLAGS) -fsyntax-only -c $*.f90
	@touch $@

$(ALL_DEPS) : $(DEPEND_DIR)/%.dep : | $(DEPEND_DIR)
	@echo MAKEDEP $*
	@$(MAKEDEP) $^ > $@

# Extra dependencies

ifneq ($(MAKECMDGOALS),clean)
  -include $(ALL_DEPS)
endif

##########
# Linking
##########

ifeq ($(shell uname -s),Darwin)
  LIB_LDFLAGS = -shared -install_name ${GYRE_DIR}/lib/$@
  EXEC_LDFLAGS = 
else			   
  LIB_LDFLAGS = -Wl,--no-warn-execstack -shared -Wl,-rpath -Wl,${GYRE_DIR}/lib
  EXEC_LDFLAGS = -Wl,--no-warn-execstack -Wl,-rpath -Wl,${GYRE_DIR}/lib
endif

LDFLAGS = `lapack_link` `hdf5_link` `lapack95_link`

ifeq ($(CRMATH),yes)
  export MESASDK_MATH_SLOT=crmath
  LDFLAGS := $(LDFLAGS) `crmath_link` `crlibm_link`
else
  export MESASDK_MATH_SLOT=default
  LDFLAGS := $(LDFLAGS) `crmath_link`
endif

# Rules

$(LIB_TARGETS) :
	@echo LD-LIB $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIB_LDFLAGS)

$(EXEC_TARGETS) : 
	@echo LD $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(EXEC_LDFLAGS)

#######################
# Secondary Expansions
#######################

.SECONDEXPANSION :

$(ALL_TARGETS) : $$(subst .fypp,.o,$$($$@_SRCS))

$(ALL_DEPS) : $(DEPEND_DIR)/%.dep : $$(subst .fypp,.f90,$$($$*_SRCS))

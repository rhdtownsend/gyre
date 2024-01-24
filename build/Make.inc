# File     : Make.inc
# Purpose  : standard makefile rules, targets

#########
# Checks
#########

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

############
# Variables
############

ifeq ($(shell uname -s),Darwin)
  LIB_SUFFIX = .dylib
else			   
  LIB_SUFFIX = .so
endif

DEPEND_DIR = .depend

BIN_FILES = $(BIN_TARGETS)
LIB_FILES = $(addsuffix $(LIB_SUFFIX),$(LIB_TARGETS))

ALL_SRCS = $(foreach target,$(BIN_TARGETS) $(LIB_TARGETS),$($(target)_SRCS) $($(target)_USES))

ALL_DPPS = $(addprefix $(DEPEND_DIR)/,$(ALL_SRCS:%.fypp=%.dpp))
ALL_DEPS = $(addprefix $(DEPEND_DIR)/,$(BIN_TARGETS:%=%.dep) $(LIB_TARGETS:%=%.dep))

all : $(BIN_FILES) $(LIB_FILES)

install : bin-install lib-install inc-install

bin-install : $(BIN_FILES)
	@for file in $^; do \
	    echo CP $${file} ${BIN_DIR}; \
	    cp $${file} ${BIN_DIR}; \
	done

lib-install : $(LIB_FILES)
	@for file in $^; do \
	    echo CP $${file} ${LIB_DIR}; \
	    cp $${file} ${LIB_DIR}; \
	done

inc-install : $(INC_FILES)
	@for file in $^; do \
	    echo CP $${file} ${INC_DIR}; \
	    cp $${file} ${INC_DIR}; \
	done

$(DEPEND_DIR) :
	mkdir -p $@

clean : almostclean
	rm -f $(DEPEND_DIR)/*

almostclean :
	rm -f $(BIN_FILES) $(LIB_FILES) *.f90 *.mod *.smod *.o 

.PHONY: all clean almostclean

.SUFFIXES:

########
# Paths
########

vpath %.fypp $(SRC_DIRS)
vpath %.inc $(SRC_DIRS)

#################
# Preproccessing
#################

# Flags & files

FPP = ./fypp_deps

FPPFLAGS := $(FPPFLAGS) \
            $(addprefix -I,$(SRC_DIRS)) \
            $(foreach opt,DEBUG OMP FPE IEEE,$(if $(filter yes,$($(opt))),-D$(opt)=1,-D$(opt)=0))

ifeq ($(DEBUG),yes)
  FPPFLAGS := $(FPPFLAGS) -DLOG_LEVEL=\'DEBUG\'
else
  FPPFLAGS := $(FPPFLAGS) -DLOG_LEVEL=\'INFO\'
endif

# Rules

%.f90 : %.fypp $(DEPEND_DIR)/%.dpp | $(DEPEND_DIR)
	@$(FPP) $(FPPFLAGS) -d $(DEPEND_DIR)/$*.dpp $< $@

$(ALL_DPPS) : 

# Include generated dependencies from .dpp files

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

# Include generated dependencies from .dep files

ifneq ($(MAKECMDGOALS),clean)
  -include $(ALL_DEPS)
endif

##########
# Linking
##########

ifeq ($(shell uname -s),Darwin)
  LIB_LDFLAGS = -shared -install_name $(LIB_DIR)/$@
  BIN_LDFLAGS = -Wl,
else			   
  LIB_LDFLAGS = -Wl,--no-warn-execstack -shared -Wl,-rpath -Wl,$(LIB_DIR)
  BIN_LDFLAGS = -Wl,--no-warn-execstack -Wl,-rpath -Wl,$(LIB_DIR)
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

$(BIN_FILES) :
	@echo LD $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(BIN_LDFLAGS)

$(LIB_FILES) :
	@echo LD-LIB $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIB_LDFLAGS)

#######################
# Secondary Expansions
#######################

.SECONDEXPANSION :

$(BIN_FILES) : % : $$(subst .fypp,.o,$$($$*_SRCS)) $$(addsuffix $(LIB_SUFFIX),$$($$*_LIBS))

$(LIB_FILES) : %$(LIB_SUFFIX) : $$(subst .fypp,.o,$$($$*_SRCS)) $$(addsuffix $(LIB_SUFFIX),$$($$*_LIBS))

$(ALL_DEPS) : $(DEPEND_DIR)/%.dep : $$(subst .fypp,.f90,$$($$*_SRCS) $$($$*_USES))

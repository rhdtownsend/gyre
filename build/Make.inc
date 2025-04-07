# File     : Make.inc
# Purpose  : standard makefile rules, targets

#########
# Checks
#########

MINIMUM_VERSION = 21.4.1

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
   STATIC_LIB_SUFFIX = .a
   SHARED_LIB_SUFFIX = .dylib
else
   STATIC_LIB_SUFFIX = .a
   SHARED_LIB_SUFFIX = .so
endif

DEPEND_DIR = .depend

BIN_FILES = $(BIN_TARGETS)

STATIC_LIB_FILES = $(addsuffix $(STATIC_LIB_SUFFIX),$(LIB_TARGETS))
SHARED_LIB_FILES = $(addsuffix $(SHARED_LIB_SUFFIX),$(LIB_TARGETS))

ifeq ($(SHARED),yes)
   LIB_FILES = $(STATIC_LIB_FILES) $(SHARED_LIB_FILES)
   LIB_SUFFIX = $(SHARED_LIB_SUFFIX)
else
   LIB_FILES = $(STATIC_LIB_FILES)
   LIB_SUFFIX = $(STATIC_LIB_SUFFIX)
endif

ALL_SRCS = $(foreach target,$(BIN_TARGETS) $(LIB_TARGETS),$($(target)_SRCS) $($(target)_USES))

ALL_DPPS = $(addprefix $(DEPEND_DIR)/,$(ALL_SRCS:%.fypp=%.dpp))
ALL_DEPS = $(addprefix $(DEPEND_DIR)/,$(BIN_TARGETS:%=%.dep) $(LIB_TARGETS:%=%.dep))

##########
# Targets
##########

install : $(addprefix $(BIN_DIR)/,$(BIN_FILES)) \
          $(addprefix $(LIB_DIR)/,$(LIB_FILES)) \
          $(addprefix $(INC_DIR)/,$(INC_FILES))
	@:

build : $(BIN_FILES) $(LIB_FILES)
	@:

clean :
	@rm -rf $(BIN_FILES) $(STATIC_LIB_FILES) $(SHARED_LIB_FILES) \
                $(INC_FILES) $(DEPEND_DIR) *.f90 *.mod *.smod *.o

$(BIN_DIR)/% : %
	@echo CP $* $@
	@cp $< $@

$(LIB_DIR)/% : %
	@echo CP $* $@
	@cp $< $@

$(INC_DIR)/% : %
	@echo CP $* $@
	@cp $< $@

$(DEPEND_DIR) :
	@mkdir -p $@

.PHONY: install build clean

.SUFFIXES:

########
# Paths
########

vpath %.fypp $(SRC_DIRS)
vpath %.inc $(SRC_DIRS)

#################
# Preproccessing
#################

# Flags

FPP = PYTHONPATH=$(MESASDK_ROOT)/lib/python ./fypp_deps

FPPOPTS = DEBUG OMP FPE IEEE

FPPFLAGS += --line-numbering $(addprefix -I,$(SRC_DIRS)) \
            $(foreach opt,$(FPPOPTS),$(if $(filter yes,$($(opt))),-D$(opt)=1,-D$(opt)=0))

ifeq ($(DEBUG),yes)
   FPPFLAGS += -DLOG_LEVEL=\'DEBUG\'
else
   FPPFLAGS += -DLOG_LEVEL=\'INFO\'
endif

# Rules

%.f90 : %.fypp | $(DEPEND_DIR)
	@$(FPP) $(FPPFLAGS) -d $(DEPEND_DIR)/$*.dpp $< $@

$(ALL_DPPS) : 

# Include generated dependencies from .dpp files

ifneq ($(MAKECMDGOALS),clean)
   -include $(ALL_DPPS)
endif

############
# Compiling
############

# Flags

FC = gfortran

MAKEDEP = awk -f makedepf08.awk

FFLAGS += -std=f2018 -finit-real=snan -fmax-errors=25 -fPIC

ifeq ($(FPE),yes)
   FFLAGS += -ffpe-trap=invalid,overflow,zero
endif

ifeq ($(DEBUG),yes)
   FFLAGS += -Og -ggdb -fcheck=all \
             -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized \
             -finline-limit=0 -fbacktrace
else
   FFLAGS += -O2 -fno-backtrace
endif

ifeq ($(OMP),yes)
   FFLAGS += -fopenmp
endif

ifeq ($(PORTABLE),yes)
   FFLAGS += -ffp-contract=off
else
   FFLAGS += -march=native
endif

# Rules

%.o : $(DEPEND_DIR)/%.anc
	@echo FC $*
	@$(FC) $(FFLAGS) -c $*.f90

$(DEPEND_DIR)/%.anc : %.f90 | $(DEPEND_DIR)
	@$(FC) $(FFLAGS) -fsyntax-only -c $*.f90
	@touch $@

$(ALL_DEPS) : | $(DEPEND_DIR)
	@echo MAKEDEP $*
	@$(MAKEDEP) $^ > $@

# Include generated dependencies from .dep files

ifneq ($(MAKECMDGOALS),clean)
   -include $(ALL_DEPS)
endif

##########
# Linking
##########

AR = ar

ifeq ($(shell uname -s),Darwin)
   SHARED_LIB_LDFLAGS = -shared -install_name $(LIB_DIR)/$@
   BIN_LDFLAGS =
   ARFLAGS = -cr
else			   
   SHARED_LIB_LDFLAGS = -shared -Wl,-rpath -Wl,$(LIB_DIR)
   ifeq ($(SHARED),yes)
      BIN_LDFLAGS = -Wl,-rpath -Wl,$(LIB_DIR)
   else
      BIN_LDFLAGS =
   endif
   ARFLAGS = cr
endif

ifeq ($(CRMATH),yes)
   export MESASDK_MATH_SLOT=crmath
   LDFLAGS += `crmath_link` `crlibm_link`
else
   export MESASDK_MATH_SLOT=default
   LDFLAGS += `crmath_link`
endif

# Rules

$(BIN_FILES) :
	@echo LD $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(BIN_LDFLAGS)

$(STATIC_LIB_FILES) :
	@echo AR-LIB $@
	@$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB_FILES) :
	@echo LD-LIB $@
	@$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(SHARED_LIB_LDFLAGS)

#######################
# Secondary Expansions
#######################

.SECONDEXPANSION :

$(BIN_FILES) : % : $$(subst .fypp,.o,$$($$*_SRCS)) $$(addsuffix $(LIB_SUFFIX),$$($$*_LIBS))

$(STATIC_LIB_FILES) : %$(STATIC_LIB_SUFFIX) : $$(subst .fypp,.o,$$($$*_SRCS))

$(SHARED_LIB_FILES) : %$(SHARED_LIB_SUFFIX) : $$(subst .fypp,.o,$$($$*_SRCS)) $$(addsuffix $(SHARED_LIB_SUFFIX),$$($$*_LIBS))

$(ALL_DEPS) : $(DEPEND_DIR)/%.dep : $$(subst .fypp,.f90,$$($$*_SRCS) $$($$*_USES))

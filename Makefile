# File     : Makefile
# Purpose  : top-level makefile

# Build frontend executables (gyre, gyre_tides, etc)
FRONTENDS ?= yes

# Build tool executables (build_poly, eval_lambda, etc)
TOOLS ?= yes

# Build additional libraries for interfacing with other codes (e.g.,
# the gyre_mesa library)
IFACES ?= no

# Build ForUM internally. If not set to "yes", then
# you must set FORUM_LIB_DIR and FORUM_INC_DIR to
# point to where the ForUM library and module files,
# respectively, are located
FORUM ?= yes

# Enable debugging (with a performance penalty)
DEBUG ?= no

# Build & link against shared libraries
SHARED ?= yes

# Enable OpenMP parallelization
OMP ?= yes

# Enable FPE checks
FPE ?= yes

# Enable correctly rounded math functions
CRMATH ?= yes

# Enable portable math (for bit-for-bit reproducibility; setting to
# no may give a performance boost)
PORTABLE ?= yes

# Use IEEE fortran features
IEEE ?= yes

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# Export options

export FRONTENDS
export TOOLS
export IFACES
export FORUM
export DEBUG
export SHARED
export OMP
export FPE
export CRMATH
export PORTABLE
export IEEE

# Enable portable math (for bit-for-bit reproducibility; setting to
# no may give a performance boost)
PORTABLE ?= yes

# Use IEEE fortran features
export IEEE ?= yes


# General make settings

SH = /bin/bash
MAKEFLAGS += --no-print-directory

# Paths

export BIN_DIR ?= $(CURDIR)/bin
export LIB_DIR ?= $(CURDIR)/lib
export INC_DIR ?= $(CURDIR)/include

export SRC_DIR := $(CURDIR)/src
export SRC_DIRS := $(addprefix $(SRC_DIR)/,ad bvp common context diff ext	\
   frontend/gyre frontend/tools grid include interp lib math matrix mode model  \
   nad output par poly rad rot sad search tar tide tnad)

ifeq ($(FORUM),yes)
   export FORUM_LIB_DIR = $(LIB_DIR)
   export FORUM_INC_DIR = $(INC_DIR)
endif

ifeq ($(CRMATH),yes)
   SRC_DIRS += $(SRC_DIR)/math/crmath
else
   SRC_DIRS += $(SRC_DIR)/math/intrinsic
endif

# Rules

install : build | $(BIN_DIR) $(LIB_DIR) $(INC_DIR)
	@$(MAKE) -C build $@

build : install-forum
	@$(MAKE) -C build $@

clean : clean-forum
	@$(MAKE) -C build $@
	@rm -rf $(BIN_DIR) $(LIB_DIR) $(INC_DIR)

test build_ref build_ref_arch :
	@$(MAKE) --no-print-directory -C test $@

check_src :
	@$(MAKE) -C build $@

ifeq ($(FORUM),yes)

   install-forum : | $(BIN_DIR) $(LIB_DIR) $(INC_DIR)
	@$(MAKE) -C $(SRC_DIR)/forum

   clean-forum :
	@$(MAKE) -C $(SRC_DIR)/forum clean

   install-forum : TESTS = no

else

   install-forum : ;
   clean-forum : ;

endif

.PHONY: install build clean test check_src build_ref build_ref_arch install-forum clean-forum

$(BIN_DIR) $(LIB_DIR) $(INC_DIR) :
	@mkdir -p $@

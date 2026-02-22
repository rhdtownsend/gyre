# File     : Makefile
# Purpose  : top-level makefile

# Build frontend executables (gyre, gyre_tides, etc)
FRONTENDS ?= yes

# Build tool executables (build_poly, eval_lambda, etc)
TOOLS ?= yes

# Build additional libraries for interfacing with other codes (e.g.,
# the gyre_mesa library)
IFACES ?= no

# Link against an external ForUM library
#
# If set to "yes", then the build system will use pkgconf to search
# for library, with a package name speficied by EXTERNAL_FORUM_PKG.
# Otherwise, the ForUM library will be built and linked internally
EXTERNAL_FORUM ?= no
EXTERNAL_FORUM_PKG ?= forum

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
# yes may incur a small performance hit)
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
export EXTERNAL_FORUM
export EXTERNAL_FORUM_PKG
export DEBUG
export SHARED
export OMP
export FPE
export CRMATH
export PORTABLE
export IEEE

# General make settings

SH = /bin/bash
MAKEFLAGS += --no-print-directory

# Paths

export BIN_DIR ?= $(CURDIR)/bin
export LIB_DIR ?= $(CURDIR)/lib
export PKG_DIR ?= $(LIB_DIR)/pkgconfig
export INC_DIR ?= $(CURDIR)/include

export SRC_DIR := $(CURDIR)/src
export SRC_DIRS := \
   $(shell find $(SRC_DIR)/eqns -type d -print) \
   $(shell find $(SRC_DIR)/frontend -type d -print) \
   $(addprefix $(SRC_DIR)/, angular bvp common context diff eqns ext \
      grid include interp lib math matrix mode model output par poly search tar tide)

ifeq ($(CRMATH),yes)
   SRC_DIRS += $(SRC_DIR)/math/crmath
else
   SRC_DIRS += $(SRC_DIR)/math/intrinsic
endif

# Rules

install : build | $(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR)
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

ifneq ($(EXTERNAL_FORUM),yes)

   install-forum : | $(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR)
	@$(MAKE) -C $(SRC_DIR)/forum

   clean-forum :
	@$(MAKE) -C $(SRC_DIR)/forum clean

   install-forum : TESTS = no

else

   install-forum : ;
   clean-forum : ;

endif

.PHONY: install build clean test check_src build_ref build_ref_arch install-forum clean-forum

$(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR) :
	@mkdir -p $@

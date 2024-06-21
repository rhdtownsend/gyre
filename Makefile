# File     : Makefile
# Purpose  : top-level makefile

# Build tool executables
TOOLS ?= yes

# Build additional interfaces (e.g., the gyre_mesa interface)
IFACES ?= yes

# Build ForUM internally. If not set to "yes", then
# you should set FORUM_LIB_DIR and FORUM_INC_DIR to
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

# Enable correctly-rounded math
CRMATH?=yes

# Enable portable math (for bit-for-bit reproducibility; setting to
# no may give a performance boost)
PORTABLE ?= yes

# Use IEEE fortran features
IEEE ?= yes

# Enable experimental features
EXPERIMENTAL ?= no

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# General make settings

export

SH = /bin/bash
MAKEFLAGS += --no-print-directory

# Paths

BIN_DIR = $(CURDIR)/bin
LIB_DIR = $(CURDIR)/lib
INC_DIR = $(CURDIR)/include
SRC_DIR = $(CURDIR)/src

ifeq ($(FORUM),yes)
   FORUM_LIB_DIR = $(LIB_DIR)
   FORUM_INC_DIR = $(INC_DIR)
endif

SRC_DIRS = $(addprefix $(SRC_DIR)/,ad bvp cad common context diff ext	\
   extern frontend/gyre frontend/tools grid include interp lib math	\
   matrix mode model nad output par poly rad rot sad search tar tide	\
   tnad)

ifeq ($(CRMATH),yes)
   SRC_DIRS += $(SRC_DIR)/math/crmath
else
   SRC_DIRS += $(SRC_DIR)/math/intrinsic
endif

# Rules

all : install-forum
	@$(MAKE) -C build

install : all | $(BIN_DIR) $(LIB_DIR) $(INC_DIR)
	@$(MAKE) -C build install

clean : clean-forum
	@$(MAKE) -C build clean
	@rm -rf $(BIN_DIR) $(LIB_DIR) $(INC_DIR)

test build_ref build_ref_arch :
	@$(MAKE) --no-print-directory -C test $@

ifeq ($(FORUM),yes)

   install-forum : | $(BIN_DIR) $(LIB_DIR) $(INC_DIR)
	@$(MAKE) -C src/forum install

   clean-forum :
	@$(MAKE) -C src/forum clean

else

   install-forum : ;
   clean-forum : ;

endif

.PHONY: all install clean test build_ref build_ref_arch install-forum clean-forum

$(BIN_DIR) $(LIB_DIR) $(INC_DIR) :
	@mkdir -p $@

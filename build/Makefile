# File     : Makefile
# Purpose  : makefile for build subdir

# Default settings

# Enable debugging (with a performance penalty)
DEBUG ?= no

# Build tool executables
TOOLS ?= yes

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

# Targets

LIB_TARGETS = 
EXEC_TARGETS = gyre

# ifeq (${TOOLS},yes)
#   EXEC_TARGETS := ${EXEC_TARGETS} ascii_to_specint synspec_to_specint ferre_to_specint c3k_to_specint \
#                   specint_to_specint specgrid_to_photgrid specint_to_specgrid \
#                   make_bb_photgrid make_passband make_photon_passband inspect_grid
#   ifneq (${FITS_LDFLAGS},)
#     EXEC_TARGETS := ${EXEC_TARGETS} goettingen_to_specint coelho_to_specint
#   endif
# endif

# Object lists & related variables
#
# libxxx_OBJS = object files required to compile libxxx
# lbbxxx_DEPS = dependency that must be added to .libxxx.d file  
#
# prog_OBJS = object files required to compile prog
# prog_LIBS = library files required to link prog
# prog_DEPS = dependencies that must be added to .prog.d file

SRCDIRS = ad bvp common context diff ext grid interp math matrix mode model nad output par rad rot sad search tar tide tnad
SRCFILES = $(notdir $(foreach dir,${SRCDIRS},$(wildcard ../src/${dir}/*_m.fypp)))
OBJFILES = $(filter-out math_%_m.o,$(patsubst %.fypp,%.o,${SRCFILES})) $(addsuffix .o,kinds_m order_m clinalg_m hgroup_m cfunc_m system_m string_m memory_m math_m)

gyre_OBJS = gyre.o ${OBJFILES}
gyre_LIBS =
gyre_DEPS = 

gyre_tides_OBJS = gyre_tides.o ${OBJFILES}
gyre_tides_LIBS =
gyre_tides_DEPS = 

# Search subdirs

GITMODS = forum
SUBDIRS = $(filter-out ${GITMODS},$(wildcard ../src/*)) $(wildcard $(patsubst %,../src/%/src,${GITMODS}))

# Include other makefiles

include Make.inc

# Additional rules

ifeq ($(CRMATH),yes)
math_m.fypp: math_crmath_m.fypp
	@ln -sf $< $@
else
math_m.fypp: math_intrinsic_m.fypp
	@ln -sf $< $@
endif

.INTERMEDIATE: math_m.fypp

.depend: math_m.fypp

# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

TARGETS=gyre poly_to_fgong poly_to_txt

ifeq ($(TIDES),yes)
  TARGETS:=${TARGETS} gyre_force gyre_response gyre_orbit
endif

ifeq ($(CONTOUR),yes)
  TARGETS:=${TARGETS} gyre_contour
endif

ODEPACK=$(shell which odepack_link > /dev/null 2>& 1 && echo "yes" || echo "no")

ifeq ($(ODEPACK),yes)
  TARGETS:=${TARGETS} build_poly
endif

# Flags

HDF5=yes
#DEBUG=yes
CRMATH=yes

LDFLAGS=`hdf5_link`

ifeq ($(ODEPACK),yes)
  LDFLAGS:=${LDFLAGS} `odepack_link`
endif

# Targets

all : 

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

include Make.inc

# Additional rules

ifeq ($(CRMATH),yes)
gyre_math.fpp: gyre_math_crmath.fpp
	@ln -sf $< $@
else
gyre_math.fpp: gyre_math_intrinsic.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: gyre_math.fpp

.depend: gyre_math.fpp

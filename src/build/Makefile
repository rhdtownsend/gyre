# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

TARGETS=gyre poly_to_fgong poly_to_txt

ifeq ($(TIDES),yes)
  TARGETS:=${TARGETS} gyre_force gyre_response gyre_orbit
endif

ODEPACK=$(shell which odepack_link > /dev/null 2>& 1 && echo "yes" || echo "no")

ifeq ($(ODEPACK),yes)
  TARGETS:=${TARGETS} build_poly
endif

# Flags

HDF5=yes
#DEBUG=yes

ifeq ($(ODEPACK),yes)
  LDFLAGS=`hdf5_link` `odepack_link`
else
  LDFLAGS=`hdf5_link`
endif

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

include Make.inc

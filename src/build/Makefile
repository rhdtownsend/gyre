# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

ODEPACK=$(shell which odepack_link > /dev/null 2>& 1 && echo "yes" || echo "no")

ifeq ($(ODEPACK),yes)
  TARGETS=gyre gyre_force poly_to_fgong build_poly poly_to_txt
else
  TARGETS=gyre gyre_force poly_to_fgong poly_to_txt
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

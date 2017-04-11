# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

ODEPACK=$(shell which madsdk_odepack_link >& /dev/null && echo "yes" || echo "no")

ifeq ($(ODEPACK),yes)
  TARGETS=gyre poly_to_fgong build_poly gyre_force
else
  TARGETS=gyre poly_to_fgong gyre_force
endif   

# Flags

HDF5=yes

ifeq ($(ODEPACK),yes)
  LDFLAGS=`mesasdk_hdf5_link` `madsdk_odepack_link`
else
  LDFLAGS=`mesasdk_hdf5_link`
endif

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../diff:../bvp:../interp:../model:../matrix:../search:../ext:../frontend:../ad:../rad:../nad:../rot:../tar:../mode:../grid:../frontend:../output:../poly

include Make.inc

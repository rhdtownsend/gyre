# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

TARGETS=gyre build_poly poly_to_fgong

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link` `madsdk_odepack_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../diff:../bep:../interp:../model:../matrix:../search:../ext:../frontend:../ad:../rad:../nad:../rot:../trad:../mode:../sol:../grid:../frontend:../output:../poly

include Make.inc

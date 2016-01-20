# File     : Makefile
# Purpose  : makefile for src/ad subdir

# Variables

TARGETS=build_poly poly_to_fgong

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link` `madsdk_odepack_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../model:../rot:../trad:../interp:../grid

include ../Make.inc

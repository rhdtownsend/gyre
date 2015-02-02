# File     : Makefile
# Purpose  : makefile for src/ad subdir

# Variables

TARGETS=build_poly_hom build_poly

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link` `madsdk_odepack_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../model:../rot:../trad

include ../Make.inc

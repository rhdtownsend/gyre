# File     : Makefile
# Purpose  : makefile for src/ad subdir

# Variables

TARGETS=gyre_ad

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../common:../ivp:../bvp:../model:../matrix:../ext:../frontend:../rad:../rad/bound:../rad/jacobian:./bound:./jacobian

include ../Make.inc

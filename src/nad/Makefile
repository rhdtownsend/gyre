# File     : Makefile
# Purpose  : makefile for src/nad subdir

# Variables

TARGETS=gyre_nad

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../common:../ivp:../bvp:../model:../matrix:../ext:../frontend:../ad:../ad/jacobian:../ad/bound:../rad:../rad/jacobian:../rad/bound:jacobian:bound

include ../Make.inc

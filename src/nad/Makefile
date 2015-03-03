# File     : Makefile
# Purpose  : makefile for src/nad subdir

# Variables

TARGETS=gyre_nad gyre_nad_map gyre_nad_contour

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../ivp:../bvp:../interp:../model:../matrix:../search:../ext:../frontend:../ad:../ad/jacobian:../ad/bound:../rad:../rad/jacobian:../rad/bound:../rot:../trad

include ../Make.inc

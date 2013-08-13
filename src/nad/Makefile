# File     : Makefile
# Purpose  : makefile for src/nad subdir

# Variables

TARGETS=gyre_nad

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=${TOPDIR}/extern/core:../common:../ivp:../bvp:../coeffs:../matrix:../frontend:../ad:../rad

include ${TOPDIR}/Make.inc

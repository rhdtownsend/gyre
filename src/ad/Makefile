# File     : Makefile
# Purpose  : makefile for src/ad subdir

# Variables

TARGETS=gyre_ad

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=${TOPDIR}/extern/core:../common:../ivp:../bvp:../coeffs:../matrix:../frontend:../rad

include ${TOPDIR}/Make.inc

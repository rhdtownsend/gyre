# File     : Makefile
# Purpose  : makefile for src/trad subdir

# Variables

TARGETS=build_trad_fit

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../common

include ../build/Make.inc

# Disable floating-point checks

FFLAGS:=${filter-out ${FFLAGS},-ffpe-trap=%}


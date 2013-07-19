# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

TARGETS=gyre_nad gyre_nad_map

# Flags

OMP=yes
MPI=no
HDF5=yes
DEBUG=yes
DOUBLE_PRECISION=yes

LDFLAGS=`mesasdk_hdf5_link` `mesasdk_lapack95_link`

FPX3FLAGS=-DGFORTRAN_PR_56052 -DGFORTRAN_PR_56872

ifdef ($DEBUG,yes)
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=DEBUG
else
  FPX3FLAGS:=${FPX3FLAGS} -DLOG_LEVEL=INFO
endif

# Include other makefiles

SUBDIRS=../../core:../common:../ivp:../bvp:../coeffs:../matrix:../ad:../frontend

include ../../build/Make.inc

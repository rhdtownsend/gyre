# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

TARGETS=gyre_nad

# Flags

OMP=true
#MPI=true
HDF5=true
DEBUG=true
DOUBLE_PRECISION=true

LDFLAGS=`mesasdk_hdf5_link` `mesasdk_lapack95_link`

FPX3FLAGS=-DGFORTRAN_PR_56052 -DGFORTRAN_PR_56872

# Include other makefiles

SUBDIRS=../../core:../common:../ivp:../bvp:../coeffs:../matrix:../ad:../frontend

include ../../build/Make.inc

# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

TARGETS=gyre_nad

# Flags

OMP=true
#MPI=true
DEBUG=true
DOUBLE_PRECISION=true

LDFLAGS=`madsdk_hdf5_link` `madsdk_lapack95_link`

# Include other makefiles

SUBDIRS=../../core:../common:../ivp:../bvp:../coeffs:../matrix:../ad

include ../../build/Make.inc

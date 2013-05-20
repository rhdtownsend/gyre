# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

TARGETS=gyre_lib

# Rules

gyre_lib.a : gyre_lib.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar cr gyre_lib.a *.o

# Flags

OMP=true
#MPI=true
DEBUG=true
DOUBLE_PRECISION=true

FPX3FLAGS=-DGFORTRAN_PR_56052 -DGFORTRAN_PR_56872

# Include other makefiles

SUBDIRS=../../core:../common:../ivp:../bvp:../coeffs:../matrix:../ad:../frontend

include ../../build/Make.inc

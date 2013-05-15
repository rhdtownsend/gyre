# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

TARGETS=gyre_mesa

# Rules

gyre_mesa.a : gyre_mesa.o
	@echo AR gyre_mesa.a
	@rm -f gyre_mesa.a
	@ar cr gyre_mesa.a *.o

# Flags

OMP=true
#MPI=true
DEBUG=true
DOUBLE_PRECISION=true

FPX3FLAGS=-DGFORTRAN_PR_56052 -DGFORTRAN_PR_56872

# Include other makefiles

SUBDIRS=../../core:../common:../ivp:../bvp:../coeffs:../matrix:../ad

include ../../build/Make.inc

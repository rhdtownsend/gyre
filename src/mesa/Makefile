# File     : Makefile
# Purpose  : makefile for src/mesa subdir

# Variables

TARGETS=gyre_lib.a test_gyre_lib

# Include other makefiles

SUBDIRS=../extern/core:../common:../ivp:../bvp:../matrix:../model:../ext:../frontend:../rad:../rad/bound:../rad/jacobian:../ad:../ad/bound:../ad/jacobian

include ../Make.inc

# Rules

gyre_lib.a : gyre_lib.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar crs gyre_lib.a *.o


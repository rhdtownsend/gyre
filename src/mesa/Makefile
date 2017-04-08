# File     : Makefile
# Purpose  : makefile for src/mesa subdir

# Variables

TARGETS=gyre_lib.a test_gyre_lib

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../common:../diff:../bvp:../interp:../model:../matrix:../search:../ext:../ad:../rad:../nad:../rot:../tar:../mode:../grid:../output

include ../build/Make.inc

# Rules

gyre_lib.a : gyre_lib.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar crs gyre_lib.a *.o


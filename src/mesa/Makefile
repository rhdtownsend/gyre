# File     : Makefile
# Purpose  : makefile for src/mesa subdir

# Variables

TARGETS=gyre_lib.a test_gyre_lib

# Flags

HDF5=yes
#DEBUG=yes
CRMATH=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../common:../math:../context:../par:../diff:../bvp:../interp:../model:../matrix:../search:../ext:../ad:../rad:../nad:../rot:../tar:../mode:../grid:../output

include ../build/Make.inc

# Rules

gyre_lib.a : gyre_lib.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar crs gyre_lib.a *.o

ifeq ($(CRMATH),yes)
gyre_math.fpp: gyre_math_crmath.fpp
	@ln -sf $< $@
else
gyre_math.fpp: gyre_math_intrinsic.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: gyre_math.fpp

.depend: gyre_math.fpp

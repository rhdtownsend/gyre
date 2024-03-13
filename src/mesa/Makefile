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

SUBDIRS=../core:../common:../math:../context:../par:../diff:../bvp:../interp:../model:../matrix:../search:../ext:../ad:../rad:../nad:../tnad:../rot:../tar:../mode:../grid:../output

include ../build/Make.inc

# Rules

gyre_lib.a : gyre_lib_m.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar crs gyre_lib.a *.o

ifeq ($(CRMATH),yes)
math_m.fpp: math_crmath_m.fpp
	@ln -sf $< $@
else
math_m.fpp: math_intrinsic_m.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: math_m.fpp

.depend: math_m.fpp

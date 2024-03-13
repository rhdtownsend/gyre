# File     : Makefile
# Purpose  : makefile for src/interp subdir

# Variables

TARGETS=test_interp

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../core:../common:../math:../matrix

include ../build/Make.inc

# Additional rules

ifeq ($(CRMATH),yes)
math_m.fpp: math_crmath_m.fpp
	@ln -sf $< $@
else
math_m.fpp: math_intrinsic_m.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: math_m.fpp

.depend: math_m.fpp

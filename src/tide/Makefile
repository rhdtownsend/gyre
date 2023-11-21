# File     : Makefile
# Purpose  : makefile for src/tide subdir

# Variables

TARGETS=eval_tidal_coeff

# Flags

HDF5=yes
CRMATH=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

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

# File     : Makefile
# Purpose  : makefile for src/trad subdir

# Variables

TARGETS=build_tar_fit map_tar_fit

# Include other makefiles

SUBDIRS=../core:../common:../math

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

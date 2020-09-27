# File     : Makefile
# Purpose  : makefile for src/trad subdir

# Variables

TARGETS=build_tar_fit map_tar_fit

# Include other makefiles

SUBDIRS=../extern/core:../common:../math

include ../build/Make.inc

# Additional rules

ifeq ($(CRMATH),yes)
gyre_math.fpp: gyre_math_crmath.fpp
	@ln -sf $< $@
else
gyre_math.fpp: gyre_math_intrinsic.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: gyre_math.fpp

.depend: gyre_math.fpp

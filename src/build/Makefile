# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

TARGETS=gyre gyre_tides poly_to_fgong poly_to_txt poly_to_parfait eval_lambda

ODEPACK=$(shell which odepack_link > /dev/null 2>& 1 && echo "yes" || echo "no")

ifeq ($(ODEPACK),yes)
  TARGETS:=${TARGETS} build_poly
endif

# Flags

ifeq ($(ODEPACK),yes)
  LDFLAGS:=${LDFLAGS} `odepack_link`
endif

# Targets

all : 

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../tnad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

include Make.inc

# Additional rules

ifeq ($(MPI),yes)
core_parallel.o : F9XFLAGS:=${F9XFLAGS} -std=legacy
endif

ifeq ($(CRMATH),yes)
gyre_math.fpp: gyre_math_crmath.fpp
	@ln -sf $< $@
else
gyre_math.fpp: gyre_math_intrinsic.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: gyre_math.fpp

.depend: gyre_math.fpp

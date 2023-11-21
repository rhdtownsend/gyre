# File     : Makefile
# Purpose  : makefile for build subdir

# Variables

TARGETS=gyre gyre_tides poly_to_fgong poly_to_txt eval_lambda

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

SUBDIRS=../core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../tnad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

include Make.inc

# Additional rules

ifeq ($(MPI),yes)
core_parallel.o : F9XFLAGS:=${F9XFLAGS} -std=legacy
endif

ifeq ($(CRMATH),yes)
math_m.fpp: math_crmath_m.fpp
	@ln -sf $< $@
else
math_m.fpp: math_intrinsic_m.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: math_m.fpp

.depend: math_m.fpp

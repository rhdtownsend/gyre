# File     : Makefile
# Purpose  : makefile for src/mesa subdir

# Variables

TARGETS=gyre_lib.a test_gyre_lib

# Flags

HDF5=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../ivp:../bvp:../interp:../model:../matrix:../search:../ext:../frontend:../rad:../rad/bound:../rad/jacobian:../rot:../trad:../ad:../nad
#SUBDIRS=../extern/core:../extern/astro:../common:../ivp:../bvp:../matrix:../model:../search:../ext:../frontend:../rad:../rad/bound:../rad/jacobian:../ad:../ad/bound:../ad/jacobian:../nad:../nad/bound:../nad/jacobian:../rot:../trad

include ../Make.inc

# Rules

gyre_lib.a : gyre_lib.o
	@echo AR gyre_lib.a
	@rm -f gyre_lib.a
	@ar crs gyre_lib.a *.o


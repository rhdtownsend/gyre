# File     : Makefile
# Purpose  : makefile for src/tide subdir

# Variables

TARGETS=gyre_orbit gyre_displace

# Flags

HDF5=yes
DEBUG=yes

LDFLAGS=`mesasdk_hdf5_link`

# Include other makefiles

SUBDIRS=../extern/core:../extern/astro:../common:../math:../par:../diff:../bvp:../interp:../model:../matrix:../search:../context:../ext:../frontend:../ad:../rad:../nad:../sad:../rot:../tar:../tide:../mode:../grid:../frontend:../output:../poly

include ../build/Make.inc


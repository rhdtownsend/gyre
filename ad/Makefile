# $Id: Makefile 121 2011-06-19 03:40:55Z townsend $

# Variables

#TARGETS=gyre_ad gyre_ad_map build_poly_hom time_gyre_ad dump_poly
TARGETS=gyre_ad build_poly_hom #time_gyre_ad dump_poly

# Flags

OMP=true
#MPI=true
DEBUG=true
DOUBLE_PRECISION=true

LDFLAGS=`mesasdk_hdf5_link` `mesasdk_lapack95_link`

# Include other makefiles

SUBDIRS=core:common:ivp:coeffs

include ../Make.inc

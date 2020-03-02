# Variables

TARGETS=test_astro_hough eval_hough

# Flags

OMP=yes
MPI=no
DEBUG=yes
HDF5=yes

DOUBLE_PRECISION=no

FPX3FLAGS:=${FPX3FLAGS} -DHDF5

LDFLAGS=`mesasdk_hdf5_link` `mesasdk_lapack95_link`

# Include other makefiles

SUBDIRS=../core

include ../build/Make.inc


# File     : Makefile
# Purpose  : top-level makefile

# Variables

TOPDIR=${CURDIR}
BINDIR=${TOPDIR}/bin

# Rules

all :
	make TOPDIR=${TOPDIR} BINDIR=${BINDIR} -w -C src

clean :
	make TOPDIR=${TOPDIR} -w -C src clean

test :
	make TOPDIR=${TOPDIR} BINDIR=${BINDIR} -w -C test

build_ref :
	make TOPDIR=${TOPDIR} BINDIR=${BINDIR} -w -C test build_ref

.PHONY: all test clean

# File     : Makefile
# Purpose  : top-level makefile

# Variables

BINDIR=${CURDIR}/bin

# Rules

all :
	@make BINDIR=${BINDIR} -w -C src install

test :
	@make BINDIR=${BINDIR} -w -C test

build_ref :
	@make BINDIR=${BINDIR} -w -C test build_ref

clean :
	@make -w -C src clean
	rm -f ${BINDIR}/*

.PHONY: all test build_ref clean

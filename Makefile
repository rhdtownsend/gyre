# File     : Makefile
# Purpose  : top-level makefile

# Variables

BINDIR=${CURDIR}/bin

# Rules

all :
	@mkdir -p ${BINDIR}
	@${MAKE} BINDIR=${BINDIR} -w -C src install

test :
	@${MAKE} BINDIR=${BINDIR} -w -C test

build_ref :
	@${MAKE} BINDIR=${BINDIR} -w -C test build_ref

clean :
	@${MAKE} -w -C src clean
	rm -f ${BINDIR}/*

.PHONY: all test build_ref clean

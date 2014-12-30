# File     : Makefile
# Purpose  : top-level makefile

# Variables

BINDIR=${CURDIR}/bin

# Rules

all :
	@mkdir -p ${BINDIR}
	@${MAKE} BINDIR=${BINDIR} -w -C src install

test :
	@${MAKE} BINDIR=${BINDIR} -w -C test $@

build_ref build_ref_arch :
	@${MAKE} BINDIR=${BINDIR} -w -C test $@

clean :
	@${MAKE} -w -C src $@
	rm -f ${BINDIR}/*

.PHONY: all test build_ref build_ref_arch clean

# File     : Makefile
# Purpose  : top-level makefile

# Variables

BINDIR=${CURDIR}/bin

# Rules

all :
	@mkdir -p ${BINDIR}
	@${MAKE} --no-print-directory BINDIR=${BINDIR} -C src install

test :
	@${MAKE} --no-print-directory BINDIR=${BINDIR} -C test $@

build_ref build_ref_arch :
	@${MAKE} BINDIR=${BINDIR} -w -C test $@

clean almostclean :
	@${MAKE} -w -C src $@
	rm -f ${BINDIR}/*

.PHONY: all test build_ref build_ref_arch clean

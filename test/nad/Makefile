# File     : Makefile
# Purpose  : makefile for test/ad subdir

# Variables

TESTDIRS=mesa/bcep mesa/rgb mesa/spb

# Rules

test :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ${TOPDIR}/test/scripts/test_gyre.sh ${BINDIR}/gyre_nad gyre_nad.in; \
            if [ $$? -ne 0 ]; then \
                exit 1; \
            fi; \
            cd -; \
        done

.PHONY: test

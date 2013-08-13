# File     : Makefile
# Purpose  : makefile for test/ad subdir

# Variables

TESTDIRS=b3/degroote hom mesa/bcep mesa/rgb mesa/spb poly/hom

# Rules

test :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ${TOPDIR}/test/scripts/test_gyre.sh ${BINDIR}/gyre_ad gyre_ad.in; \
            if [ $$? -ne 0 ]; then \
                exit 1; \
            fi; \
            cd -; \
        done

.PHONY: test

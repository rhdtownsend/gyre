# File     : Makefile
# Purpose  : makefile for test/ad subdir

# Variables

TESTDIRS=hom poly/hom mesa/bcep mesa/rgb mesa/spb 

# Rules

test :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ${TOPDIR}/test/scripts/test_gyre.sh ${BINDIR}/gyre_ad ${CURDIR}/$$DIR/gyre_ad.in "`cat label`"; \
            if [ $$? -ne 0 ]; then \
                exit 1; \
            fi; \
            cd ${CURDIR}; \
        done

.PHONY: test

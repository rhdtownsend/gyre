# File     : Makefile
# Purpose  : makefile for test/ad subdir

# Variables

TESTDIRS=hom poly/hom omp doppler mesa/bcep mesa/rgb mesa/spb solar

# Rules

test :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./test.sh; \
            cd ${CURDIR}; \
	done

build_ref :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./gyre_ad gyre_ad.in; \
            for file in *.txt; do mv $$file ref/$$file; done; \
            cd ${CURDIR}; \
        done

build_ref_arch :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./gyre_ad gyre_ad.in; \
            for file in *.txt; do mv $$file ref/$$file.`uname -m`-`uname -s`; done; \
            cd ${CURDIR}; \
        done

.PHONY: test build_ref build_ref_arch

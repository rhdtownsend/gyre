# File     : Makefile
# Purpose  : makefile for test/nad subdir

# Variables

TESTDIRS=mesa/bcep mesa/rgb mesa/spb

# Rules

test :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./test.sh; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
            cd ${CURDIR}; \
        done

build_ref :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./gyre_nad gyre_nad.in; \
            for file in *.txt; do mv $$file ref/$$file; done; \
            cd ${CURDIR}; \
        done

build_ref_arch :
	@for DIR in ${TESTDIRS}; do \
            cd $$DIR; \
            ./gyre_nad gyre_nad.in; \
            for file in *.txt; do mv $$file ref/$$file.`uname -m`-`uname -s`; done; \
            cd ${CURDIR}; \
        done

.PHONY: test build_ref

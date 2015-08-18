# File     : Makefile
# Purpose  : makefile for num subdir

# Variables

SUBDIRS=omp band doppler trad cowling vars-jcd vars-mix

# Rules

test :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR}; \
        done

build_ref build_ref_arch :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR} $@; \
        done

.PHONY: test build_ref build_ref_arch

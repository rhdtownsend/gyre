# File     : Makefile
# Purpose  : makefile for poly subdir

# Variables

SUBDIRS=0.0 0.0+0.0 1.5 1.5-no-ctr 3.0 3.0-g 3.0+3.0 4.0

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

# File     : Makefile
# Purpose  : makefile for mesa subdir

# Variables

SUBDIRS=bcep spb spb-rossby rgb dcep solar

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

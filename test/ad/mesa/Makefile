# File     : Makefile
# Purpose  : makefile for mesa subdir

# Variables

SUBDIRS=bcep spb spb-qad spb-linear-interp rgb rgb-isolate-gamma rgb-isolate-pi dcep

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

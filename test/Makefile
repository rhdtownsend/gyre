# File     : Makefile
# Purpose  : makefile for test subdir

# Variables

SUBDIRS=ad nad tides

# Rules

ifneq ($(PORTABLE),yes)
   $(warning Warning: PORTABLE is not set to "yes", some tests may fail)
endif

test :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR}; \
        done

build_ref build_ref_arch :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR} $@; \
        done

.PHONY: test build_ref build_ref_arch

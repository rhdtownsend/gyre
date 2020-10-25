# File     : Makefile
# Purpose  : makefile for test subdir

# Variables

SUBDIRS=ad nad

ifeq ($(TIDE),yes)
  SUBDIRS:=${SUBDIRS} tide
endif

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

# File     : Makefile
# Purpose  : makefile for test subdir

# Variables

SUBDIRS=ad nad

# Rules

test :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR}; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

build_ref build_ref_arch :
	@for DIR in ${SUBDIRS}; do \
            make -C $${DIR} $@; \
        done

.PHONY: test build_ref

# File     : Makefile
# Purpose  : makefile for test subdir

# Variables

SUBDIRS=ad nad

# Rules

test :
	@for DIR in ${SUBDIRS}; do \
            make BINDIR=${BINDIR} -w -C $${DIR}; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

build_ref :
	@for DIR in ${SUBDIRS}; do \
            make BINDIR=${BINDIR} -w -C $${DIR} build_ref; \
        done

.PHONY: test build_ref

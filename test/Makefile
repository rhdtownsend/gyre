# File     : Makefile
# Purpose  : makefile for test subdir

# Variables

SUBDIRS=ad nad

# Rules

test :
	@for DIR in ${SUBDIRS}; do \
            make TOPDIR=${TOPDIR} BINDIR=${BINDIR} -w -C $${DIR}; \
        done

.PHONY: test

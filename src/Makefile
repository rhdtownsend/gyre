# File     : Makefile
# Purpose  : makefile for src subdir

# Variables

SUBDIRS=ad nad

# Rules

all :
	@for DIR in ${SUBDIRS}; do \
            make TOPDIR=${TOPDIR} BINDIR=${BINDIR} -w -C $${DIR}; \
        done

clean :
	@for DIR in ${SUBDIRS}; do \
            make TOPDIR=${TOPDIR} -w -C $${DIR} clean; \
        done

.PHONY: all clean

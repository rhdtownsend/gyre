# File     : Makefile
# Purpose  : makefile for src subdir

# Variables

SUBDIRS=ad nad

# Rules

all :
	@for DIR in ${SUBDIRS}; do \
            make -w -C $${DIR}; \
        done

install :
	@for DIR in ${SUBDIRS}; do \
            make BINDIR=${BINDIR} -w -C $${DIR} install; \
        done

clean :
	@for DIR in ${SUBDIRS}; do \
            make -w -C $${DIR} clean; \
        done

.PHONY: all clean

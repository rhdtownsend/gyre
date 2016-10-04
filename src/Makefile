# File     : Makefile
# Purpose  : makefile for src subdir

# Variables

SUBDIRS=build

# Rules

all :
	@for DIR in ${SUBDIRS}; do \
            ${MAKE} -w -C $${DIR}; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

install :
	@for DIR in ${SUBDIRS}; do \
            ${MAKE} BINDIR=${BINDIR} -w -C $${DIR} install; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

clean almostclean :
	@for DIR in ${SUBDIRS}; do \
            ${MAKE} -w -C $${DIR} $@ ; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

.PHONY: all clean almostclean

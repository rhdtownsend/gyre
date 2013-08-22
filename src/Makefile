# File     : Makefile
# Purpose  : makefile for src subdir

# Variables

SUBDIRS=ad nad

# Rules

all :
	@for DIR in ${SUBDIRS}; do \
            make -w -C $${DIR}; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

install :
	@for DIR in ${SUBDIRS}; do \
            make BINDIR=${BINDIR} -w -C $${DIR} install; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

clean :
	@for DIR in ${SUBDIRS}; do \
            make -w -C $${DIR} clean; \
            if [ $$? -ne 0 ]; then \
                exit $$?; \
            fi; \
        done

.PHONY: all clean

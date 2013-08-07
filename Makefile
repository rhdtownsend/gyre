# File     : Makefile
# Purpose  : top-level makefile

# Variables

MAKEDIRS=ad nad

TOPDIR=${CURDIR}
BINDIR=${TOPDIR}/bin

# Rules

all :
	@for DIR in $(addprefix src/,${MAKEDIRS}); do \
            make BINDIR=$${BINDIR} -w -C $${DIR}; \
        done

test :
	@for DIR in $(addprefix test/,${MAKEDIRS}); do \
            make BINDIR=$${BINDIR} -w -C $${DIR}; \
        done

clean :
	@for DIR in $(addprefix src/,${MAKEDIRS}); do \
            make -w -C $${DIR} clean; \
        done
	@rm -f ${BINDIR}/*

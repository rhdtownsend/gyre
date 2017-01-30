# File     : Makefile
# Purpose  : makefile for num subdir

# Variables

SUBDIRS=omp band delta-p delta-g doppler tar cowling vars-dziem vars-jcd vars-mix vars-lagp hdf

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

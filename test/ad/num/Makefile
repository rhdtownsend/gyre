# File     : Makefile
# Purpose  : makefile for num subdir

# Variables

SUBDIRS=omp band delta cowling vars-dziem vars-jcd vars-mix vars-lagp \
        diff-colloc diff-magnus diff-mirk diff-trapz \
        lambda-sph lambda-tar-gravity lambda-tar-rossby lambda-adhoc

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

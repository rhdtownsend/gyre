# File     : Makefile
# Purpose  : makefile for analytic poly subdir

# Rules

test :
	@./test.sh

build_ref :
	./gyre gyre.in; \
        for file in summary*.h5; do \
	    if h5diff -q $$file ref/$$file; then \
                rm $$file; \
            else \
	        mv $$file ref/$$file; \
            fi; \
	done

.PHONY: test build_ref build_ref_arch

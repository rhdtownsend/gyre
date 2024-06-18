# File     : Makefile
# Purpose  : top-level makefile

# Variables

export BIN_DIR = $(CURDIR)/bin
export LIB_DIR = $(CURDIR)/lib
export INC_DIR = $(CURDIR)/include

# Rules

all : | $(BIN_DIR) $(LIB_DIR) $(INC_DIR)
	@$(MAKE) --no-print-directory -C src/forum/build
	@$(MAKE) --no-print-directory -C src/forum/build install
	@$(MAKE) --no-print-directory -C build
	@$(MAKE) --no-print-directory -C build install

test :
	@$(MAKE) --no-print-directory BINDIR=$(BINDIR) -C test $@

build_ref build_ref_arch :
	@$(MAKE) -w -C test $@

clean almostclean :
	@$(MAKE) --no-print-directory -C src/forum/build $@
	@$(MAKE) --no-print-directory -C build $@
	rm -f $(BIN_DIR)/* $(LIB_DIR)/* $(INC_DIR)/*

.PHONY: all test build_ref build_ref_arch clean almostclean

$(BIN_DIR) $(LIB_DIR) $(INC_DIR) :
	@mkdir -p $@

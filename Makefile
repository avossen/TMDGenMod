# Main make file to enter other directories and use those makefiles

include Makefile.in

SUBDIRS = $(wildcard *.d)

all: lib bin install

.PHONY: obj bin lib clean $(SUBDIRS)

# make in each subdir

obj: $(SUBDIRS)

$(SUBDIRS):
	@make -C $@ obj

lib: obj
	@make -C lib

bin: lib
	@make -C Exec.d

install: bin
	@echo
	@echo Executable TMDGen is now available in $(BASE_DIR)/bin
	@echo Libraries are now available in $(BASE_DIR)/lib
	@echo
	@echo Enjoy!
	@echo



# clean all sub directories

clean:
	@for dir in $(SUBDIRS) lib; do \
		make clean -C $$dir; \
	done

###
### ------------------------------------------------------------
### Copyright (c) 2013  S. Gliske and the University of Michigan
### This software is the proprietary information of University of Michigan.
### All Rights Reserved.
### ------------------------------------------------------------
###
### SVN revision information:
### @version $Revision$:
### @author  $Author$:
### @date    $Date$:
###

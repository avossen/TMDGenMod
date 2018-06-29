include ../Makefile.in

#GRV_DIR=$(subst DF_Wrapper,GRV,$(PWD))/grids
#CTEQ_DIR=$(subst DF_Wrapper,CTEQ6,$(PWD))
GRV_DIR=$(PWD)/GRV.d/grids
CTEQ_DIR=$(PWD)/CTEQ6.d

.cpp.o:
	$(CPP) -c $(CPP_FLAGS)  $(INCLUDE_PATH) -o $@ $*.cpp $(CFLAGS_LHAPDF) -DGRV_GRID_DIR=\"$(GRV_DIR)\" -DCTEQ_GRID_DIR=\"$(CTEQ_DIR)\"

.F.o:
	$(FORTRAN) -c $(FFLAGS)  $(INCLUDE_PATH) -o $@ $*.F -I../CTEQ6.d 

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z_].cpp)) gmc_trans_setctq6.o

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:

clean:
	@rm -f *.a *.o *.so


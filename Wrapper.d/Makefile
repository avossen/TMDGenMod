include ../Makefile.in

.cpp.o:
	$(CPP) -c $(CPP_FLAGS)  $(INCLUDE_PATH) -o $@ $*.cpp $(ROOT_CPP_FLAGS)

.F.o:
	$(FORTRAN) -c $(FFLAGS)  $(INCLUDE_PATH) -o $@ $*.F

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z_].cpp))

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:    TestWrapper

clean:
	@rm -f *.a *.o *.so TestWrapper

TestWrapper: TestWrapper.o TMDGen_Wrapper.o
	$(FORTRAN) -o TestWrapper TestWrapper.o TMDGen_Wrapper.o $(LDFLAGS) -L../lib -lTMDGen $(GSL_LIBS) $(ROOT_LD_FLAGS) -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32

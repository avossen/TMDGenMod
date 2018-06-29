include ../Makefile.in

.cpp.o:
	$(CPP) -c $(CPP_FLAGS)  $(INCLUDE_PATH) -o $@ $*.cpp

.F.o:
	$(FORTRAN) -c $(FFLAGS)  $(INCLUDE_PATH) -o $@ $*.F

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z].cpp)) fint.o

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:

clean:
	@rm -f *.a *.o *.so


include ../Makefile.in

INCLUDE_PATH += -I.

.cpp.o:
	$(CPP) -c $(CPP_FLAGS)  $(INCLUDE_PATH) -o $@ $*.cpp

.F.o:
	$(FORTRAN) -c $(FFLAGS)  $(INCLUDE_PATH) -o $@ $*.F

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z].cpp)) $(patsubst %.F,%.o,$(wildcard *.F)) 

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:

test: test.o
	$(CPP) $(LDFLAGS) $(LIBS) $(OBJECT_FILES) -o test


clean:
	@rm -f *.a *.o *.so


include ../Makefile.in

.cpp.o:
	$(CPP) -c $(CPP_FLAGS)  $(INCLUDE_PATH) -o $@ $*.cpp

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z_].cpp))

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:

clean:
	@rm -f *.a *.o *.so


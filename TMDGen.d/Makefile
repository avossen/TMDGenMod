include ../Makefile.in

INCLUDE_PATH += $(DAD_PATH)

.cpp.o:
	$(CPP) -c -o $@ $*.cpp $(CPP_FLAGS) $(ROOT_CPP_FLAGS) $(INCLUDE_PATH)  $(ROOT_CPP_FLAGS)

.PHONY: obj OBJECT_FILES

OBJECT_FILES = $(patsubst %.cpp,%.o,$(wildcard *[a-zA-Z].cpp))

all:	obj

obj:	$(OBJECT_FILES)

libs:

bin:

clean:
	@rm -f *.a *.o *.so


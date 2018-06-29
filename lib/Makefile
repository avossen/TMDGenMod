include ../Makefile.in

THESE_LIBS =  libTMDGen.a  libTMDGen_.a 

# libTMDGen.so libTMDGen_.so

OBJECTS = $(wildcard ../[a-zA-VX-Z]*.d/*.o)
REDUCED_OBJECTS = $(filter-out ../Common.d/ParseInput.o, $(OBJECTS))

all: lib

obj:
	@make -C .. obj

lib: $(THESE_LIBS)

libTMDGen.a: obj
	ar rcs libTMDGen.a  $(AR_LIBS) $(OBJECTS)

libTMDGen_.a: obj
	ar rcs libTMDGen_.a  $(AR_LIBS) $(REDUCED_OBJECTS)

libTMDGen.so: obj
	$(CPP) -shared -Wl,-soname,$@ -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS) -fPIC


libTMDGen_.so: obj
	$(CPP) -shared -Wl,-soname,$@ -o $@ $(REDUCED_OBJECTS) $(LDFLAGS) $(LIBS) -fPIC

clean:
	@rm -f *.so *.a

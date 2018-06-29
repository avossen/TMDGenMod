include ../Makefile.in

OTHER_LD_FLAGS = -L/home/sgliske/hermes/projects/cpp_code/lib -lPWC_Module -lsgUtil $(ROOT_LD_FLAGS) 
OTHER_CPP_FLAGS = -I/home/sgliske/hermes/projects/cpp_code/PWC_Module.d -I/home/sgliske/hermes/projects/cpp_code/Util.d $(ROOT_CPP_FLAGS)


OBJECTS = $(wildcard ../[a-zA-VX-Z]*.d/*.o)


.PHONY: bin EXEC_FILES

#EXEC_FILES = $(patsubst %.cpp,%,$(wildcard *[a-zA-Z].cpp))
EXEC_FILES = TMDGen_exec


OBJ_FILES = $(wildcard ../*.d/*.o)

all:	bin

bin:	clean $(EXEC_FILES) install

clean:
	@rm -f *.a *.o *.so $(EXEC_FILES)

.cpp: $(OBJECTS)
	$(CPP) $(CPP_FLAGS) -o $@ $*.cpp -I.. $(INCLUDE_PATH) $(LD_FLAGS) $(LIBS) $(GSL_LIBS) $(OBJECTS)

.PHONY: install

TMDGen_exec: TMDGen_exec.cpp ../lib/libTMDGen.a
	$(CPP) -o TMDGen TMDGen_exec.cpp $(CPP_FLAGS) -I.. $(INCLUDE_PATH) -L../lib -lTMDGen $(LD_FLAGS) $(LIBS) $(GSL_LIBS) -L../lib -lTMDGen $(ROOT_LD_FLAGS) $(ROOT_CPP_FLAGS) $(DAD_LD_FLAGS) $(LD_LHAPDF)

TMDGen_exec.cpp:
	@

Opt_D1_Spec_Ia: ../lib/libTMDGen.a
	$(CPP) $(CPP_FLAGS) -o Opt_D1_Spec_Ia Opt_D1_Spec_Ia.cpp -I.. $(INCLUDE_PATH) -L../lib -lTMDGen $(LD_FLAGS) $(LIBS) $(GSL_LIBS) $(ROOT_LD_FLAGS) $(ROOT_CPP_FLAGS) $(DAD_LD_FLAGS)

install:
	@cp TMDGen ../bin

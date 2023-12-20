SRC += $(shell find ./src -name '*.cpp')
#./PETScSolver/src 

OBJECTS = $(shell find ./src -name '*.o')
#./PETScSolver/src 

LIB_OBJECT = libpetscsolver.a
LIB_DIR = /public1/home/sch10084/petsc_solver/lib


FIM_SRC_V3 += $(SRC)

#IMPES_SRC_V3 +=$(SRC) 
FIM_OBJ_V3 = $(patsubst %.cpp, %.o, $(FIM_SRC_V3))

#IMPES_OBJ_V3 = $(patsubst %.cpp, %.o, $(IMPES_SRC_V3))

include make.inc

CLFLAGS += -I./include 
#CLFLAGS += -I./home/spring/zlj/simfast_mpi/petsc-3.6.3/include
#-Wno-c++11-compat-deprecated-writable-strings


LDFLAGS  = $(LINK_OPT)

#all: PAR_FIM
all: FIM  $(LIB_OBJECT) 
#all: FIM IMPES

FIM: $(FIM_SRC_V3) $(FIM_PROG_V3)


#IMPES: $(IMPES_SRC_V3) $(IMPES_PROG_V3)

$(FIM_PROG_V3): $(FIM_OBJ_V3)
	$(LINKER) $(FIM_OBJ_V3) -o $@ $(LDFLAGS)
	@echo 'Built executable file $@'

$(LIB_OBJECT): ${OBJECTS}
	ar rc $(LIB_OBJECT) ${OBJECTS}
	@echo 'Built lib file $@'
	mv $(LIB_OBJECT) $(LIB_DIR)

#$(IMPES_PROG_V3): $(IMPES_OBJ_V3)
#	@$(LINKER) $(IMPES_OBJ_V3) -o $@ $(LDFLAGS) 
#	@echo 'Built executable file $@'

%.o: %.cpp
	@$(CXX) $(CLFLAGS) -c -o $@ $^ 
	@echo 'Built object file $@'	

clean:
	@find . -name '*.o' -exec rm {} \;
	@rm AdjMatrixGenerator PennSimV3FIM PennSimV3WithAdjFIM ParPennSimV3FIM
#	@find . -name '.*' -exec rm {} \;

distclean:
	@find . -name '*.o' -exec rm {} \;
	@rm -f $(FIM_PROG_V3) $(IMPES_PROG_V3) 

debug:
	@echo $(CLFLAGS)

.PHONY: clean distclean debug

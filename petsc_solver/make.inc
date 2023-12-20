
PETSC_DIR = /public1/home/sch10084/petsc
PETSC_ARCH = arch-linux-c-debug
include $(PETSC_DIR)/lib/petsc/conf/variables

CC=mpiicc
#CXX=mpicxx
CXX=mpiicpc
FC=mpiifort

LINKER=$(FC)
PETSC_INCLUDE = $(CCPPFLAGS) $(PCC_FLAGS) $(CFLAGS)

 PETSC_LINK =  $(PETSC_LIB)

#PETSC_LINK =  -Wl,-rpath,/HOME/sysu_xys_3/WORKSPACE/gwench/ParPennSim-test-4 -Wl,-rpath,/HOME/sysu_xys_3/WORKSPACE/petsc-3.4.5/arch-linux2-c-opt/lib -L/WORK/sysu_xys_3/petsc-3.4.5/arch-linux2-c-opt/lib  -lpetsc -Wl,-rpath,/WORK/sysu_xys_3/petsc-3.4.5/arch-linux2-c-opt/lib -lHYPRE -Wl,-rpath,/usr/local/mpi3-dynamic/lib -L/usr/local/mpi3-dynamic/lib -Wl,-rpath,/opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64 -L/opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64 -Wl,-rpath,/opt/intel/composer_xe_2013_sp1.2.144/ipp/lib/intel64 -L/opt/intel/composer_xe_2013_sp1.2.144/ipp/lib/intel64 -Wl,-rpath,/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64 -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64 -Wl,-rpath,/opt/intel/composer_xe_2013_sp1.2.144/tbb/lib/intel64/gcc4.4 -L/opt/intel/composer_xe_2013_sp1.2.144/tbb/lib/intel64/gcc4.4 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -lmpicxx -Wl,-rpath,/HOME/sysu_xys_3/WORKSPACE/lib -L/HOME/sysu_xys_3/WORKSPACE/lib -ldl  -lmkl_intel_lp64 -lpthread -lm -Wl,-Bstatic -lmkl_sequential -lmkl_core -Wl,-Bdynamic -lX11 -lmpifort -lifport -lifcore -lmpicxx  -lmpi -limf -lsvml -lirng -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s 

##############################################################################################

#HYPRE_FLAG = -D_HYPRE_PARALLEL_SOLVER_ -DHAVE_CONFIG_H -DHYPRE_TIMING

#HYPRE_DIR = /HOME/nsfc2015_567/NSFC/gwench/soft/hypre-2.10.1
#HYPRE_INCLUDE = -I$(HYPRE_DIR)/src/hypre/include
#HYPRE_LINK = -L$(HYPRE_DIR)/src/hypre/lib -lHYPRE

##############################################################################################

LAPACKE_INCLUDE = -I//home/spring/zlj/simfast_mpi/lapack-3.6.0/LAPACKE/include/
LAPACKE_LINK = -L//home/spring/zlj/simfast_mpi/lapack-3.6.0/ -llapacke

# Set the path to your fasp and fasp4blkoil
#FASPDIR = $(HOME)/WORKSPACE/faspsolver
#FASP4BLKOILDIR= $(HOME)/WORKSPACE/fasp4blkoil

#FASP_INCLUDE = -I$(FASPDIR)/base/include -I$(FASP4BLKOILDIR)/include/ -I./FASPSolver/include/
#FASP_LIB = -L$(FASPDIR)/lib -L$(FASP4BLKOILDIR)/lib -lfasp -lfasp4blkoil

# Setup UMFPACK link options if you have compiled FASP with UMFPACK
#UMFPACKLIB_DIR= -L/usr/local/lib
#UMFPACKLIB= $(UMFPACKLIB_DIR) -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd -lsuitesparseconfig 

##############################################################################################
#Set the partition include path and library path
PARTITION_DIR=/home/spring/zlj/simfast_mpi/Partition
PARMETIS_DIR=/home/spring/zlj/simfast_mpi/parmetis-4.0.3
PARTITION_INCLUDE=-I$(PARTITION_DIR)/include -I$(PARMETIS_DIR)/include -I$(PARMETIS_DIR)/programs -I$(PARMETIS_DIR)/metis/GKlib/ -I$(PARMETIS_DIR)/parmetis -I$(PARMETIS_DIR)/metis/include
PARTITION_LINK=-lpart -lGKlib -lparmetis -lmetis -L$(PARTITION_DIR)/lib -L$(PARMETIS_DIR)/build/Linux-x86_64/libparmetis/ -L$(PARMETIS_DIR)/build/Linux-x86_64/libmetis/ -L$(PARMETIS_DIR)/metis/GKlib/build/Linux-x86_64

#############################################################################################
# Choose your BLAS and LAPACK libraries

# BLAS_LAPACK_INC = -I/usr/local/include

# Generic Linux BLAS and LAPACK setting
# BLAS_LAPACK = -lblas -llapack

# Uncomment the following line if you wish to call mkl_blas and mkl_lapack libs
# BLAS_LAPACK = -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_sequential -lmkl_core -liomp5 -lpthread

# Or uncomment the following if you are using the intel MPI compiler
#BLAS_LAPACK = -mkl

# Uncomment the following if you wish to call the system BLAS and LAPACK of Mac OS X
#BLAS_LAPACK= -framework Accelerate

#BLAS_LAPACK_FLAG = $(BLAS_LAPACK_INC) $(BLAS_LAPACK_MACRO)

##############################################################################################
# Choose c++ standard libs: -lstdc++ for linux and -lc++ for Mac OS X
STD_LIB = -lstdc++

##############################################################################################
# CLFLAGS = -g -D__SOLVER_FASP__  $(FASP_INCLUDE) $(BLAS_LAPACK_FLAG) $(PARTITION_INCLUDE) $(HYPRE_INCLUDE) $(PETSC_INCLUDE) $(HYPRE_FLAG) $(LAPACKE_INCLUDE)

# YOU CAN USE PARALLEL_DEBUG_DISPLAY MACRO TO SEE SOME VARIABLES AND -D_PETSC_VERSION_3_DOT_4_DOT_5_ FOR PETSC-3.4.5
CLFLAGS = -O2 $(BLAS_LAPACK_FLAG) $(PARTITION_INCLUDE) $(HYPRE_INCLUDE) $(PETSC_INCLUDE) $(HYPRE_FLAG) $(LAPACKE_INCLUDE)

# add option -nofor-main if your fortran compiler is ifort
LINK_OPT = -nofor-main  $(HYPRE_LINK) $(FASP_LIB) $(UMFPACKLIB) $(BLAS_LAPACK) $(STD_LIB) $(PARTITION_LINK) $(PETSC_LINK) $(LAPACKE_LINK)
#LINK_OPT =  $(HYPRE_LINK) $(FASP_LIB) $(UMFPACKLIB) $(BLAS_LAPACK) $(STD_LIB) $(PARTITION_LINK) $(PETSC_LINK) $(LAPACKE_LINK)

#SRC += $(shell find ./FASPSolver -name '*.cpp')
FIM_SRC_V3 += ./main/PennSimV3FIMmain.cpp
PAR_FIM_SRC_V3 += ./main/ParPennSimV3FIMmain.cpp
ADJ_OUTPUT_SRC_V3 += ./main/AdjGenerator.cpp
FIM_WITH_INPUT_SRC_V3 += ./main/PennSimV3FIMmainWithAdjGraph.cpp
#IMPES_SRC_V3 += ./main/PennSimV3IMPESmain.cpp
FIM_PROG_V3 = PennSimV3FIM
PAR_FIM_PROG_V3 = ParPennSimV3FIM
ADJ_OUTPUT_V3 = AdjMatrixGenerator
FIM_WITH_INPUT_PROG_V3 = PennSimV3WithAdjFIM
#IMPES_PROG_V3 = PennSimV3IMPES

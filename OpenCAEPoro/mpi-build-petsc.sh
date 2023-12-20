#!/bin/bash


# source /es01/paratera/parasoft/module.sh  
# source /es01/paratera/parasoft/oneAPI/2022.1/setvars.sh
# module load cmake/3.17.1 gcc/7.3.0-para 

export CC=mpicc
export CXX=mpicxx
#export CFLAGS += "-O3 -Xcompiler -fopenacc"
#export CXXFLAGS += "-O3 -Xcompiler -fopenacc"


# users specific directory paths
export PARMETIS_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3
export PARMETIS_BUILD_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/build/Linux-x86_64
export METIS_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/metis
export METIS_BUILD_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/build/Linux-x86_64
export PETSC_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3
export PETSC_ARCH=petsc_install
export PETSCSOLVER_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc_solver


export CPATH=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/include/:$CPATH
export CPATH=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/petsc_install/include/:/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/metis/include:/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/include:$CPATH
export CPATH=/home/gh0s7/project/OpenCAEPoro_ASC2024/lapack-3.11/CBLAS/include/:$CPATH



# install
rm -fr build; mkdir build; cd build;

echo "cmake -DUSE_PETSCSOLVER=ON -DUSE_PARMETIS=ON -DUSE_METIS=ON -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release .."
cmake -DUSE_PETSCSOLVER=ON -DUSE_PARMETIS=ON -DUSE_METIS=ON -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release ..

echo "make -j 32"
make -j 32

echo "make install"
make install

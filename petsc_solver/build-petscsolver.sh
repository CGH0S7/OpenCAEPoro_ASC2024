#!/bin/bash

# source /es01/paratera/parasoft/module.sh  
# source /es01/paratera/parasoft/oneAPI/2022.1/setvars.sh
# module load cmake/3.17.1 gcc/7.3.0-para 

export CC=mpicc
export CXX=mpicxx
export CFLAGS="-O3 -fopenacc"
export CXXFLAGS="-O3 -fopenacc"

export CPATH=/home/gh0s7/project/OpenCAEPoro_ASC2024/lapack-3.11/CBLAS/include:/home/gh0s7/project/OpenCAEPoro_ASC2024/lapack-3.11/LAPACKE/include:$CPATH
export LD_LIBRARY_PATH=/home/gh0s7/project/OpenCAEPoro_ASC2024/lapack-3.11:$LD_LIBRARY_PATH

rm -rf build 
mkdir build
cd build
cmake ..
make
mv libpetsc_solver.a ../lib/



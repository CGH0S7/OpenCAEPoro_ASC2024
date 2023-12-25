#!/bin/bash


# source /es01/paratera/parasoft/module.sh  
# source /es01/paratera/parasoft/oneAPI/2022.1/setvars.sh
# module load cmake/3.17.1 gcc/7.3.0-para 


cd src/
make clean

# configure
./configure --prefix="/home/gh0s7/project/OpenCAEPoro_ASC2024/hypre-2.28.0/install" --with-MPI --enable-shared

# compile
make -j 16
# install
make install

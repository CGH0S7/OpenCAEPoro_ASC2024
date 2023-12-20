#!/bin/bash



export PETSC_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3
export PETSC_ARCH=petsc_install

./configure  \
    --with-fortran-bindings=0 \
	--with-hypre-dir=/home/gh0s7/project/OpenCAEPoro_ASC2024/hypre-2.28.0/install \
    --with-debugging=0 \
    --with-mpi-dir=/usr/
    COPTFLAGS="-O3" \
    CXXOPTFLAGS="-O3" \


make -j 20 PETSC_DIR=/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3 PETSC_ARCH=petsc_install all
make all check    

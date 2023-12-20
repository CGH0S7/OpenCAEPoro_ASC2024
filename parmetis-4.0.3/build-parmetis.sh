#!/bin/bash

# source /es01/paratera/parasoft/module.sh  
# source /es01/paratera/parasoft/oneAPI/2022.1/setvars.sh
# module load cmake/3.17.1 gcc/7.3.0-para 


make config cc=mpicc prefix=/home/gh0s7/project/OpenCAEPoro_ASC2024/parmetis-4.0.3/parmetis-install
make -j 16
make install

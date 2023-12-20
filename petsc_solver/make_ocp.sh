#!/bin/sh

make
mv libpetsc_solver.a ../lib/
cd /public1/home/sch10084/OpenCAEPoroX/build

make -j 16


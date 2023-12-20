#!/usr/bin/python3
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--with-debugging=0',
    '--with-fortran-bindings=0',
    '--with-hypre-dir=/home/gh0s7/project/OpenCAEPoro_ASC2024/hypre-2.28.0/install',
    '--with-mpi-dir=/usr/',
    'PETSC_ARCH=petsc_install',
  ]
  configure.petsc_configure(configure_options)

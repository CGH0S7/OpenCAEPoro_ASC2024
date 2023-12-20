set(CMAKE_C_COMPILER "/public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "20.2.1.20211109")
set(CMAKE_C_PLATFORM_ID "Linux")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()




set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;dl;rt;pthread;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;irc;svml;c;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/public4/soft/oneAPI/2022.1/mpi/2021.5.0/lib/release;/public4/soft/oneAPI/2022.1/mpi/2021.5.0/lib;/opt/slurm/slurm/lib;/public4/soft/oneAPI/2022.1/mpi/2021.5.0/libfabric/lib;/public4/soft/oneAPI/2022.1/ccl/latest/lib/cpu;/public4/soft/oneAPI/2022.1/dal/latest/lib/intel64;/public4/soft/oneAPI/2022.1/tbb/latest/lib/intel64/gcc4.8;/public4/soft/oneAPI/2022.1/mkl/latest/lib/intel64;/public4/soft/oneAPI/2022.1/compiler/latest/linux/compiler/lib/intel64;/public4/soft/oneAPI/2022.1/compiler/latest/linux/lib;/public4/soft/oneAPI/2022.1/ipp/latest/lib;/public4/soft/oneAPI/2022.1/compiler/2022.0.0/linux/compiler/lib/intel64_lin;/public4/soft/gcc/9.1.0/lib/gcc/x86_64-pc-linux-gnu/9.1.0;/public4/soft/gcc/9.1.0/lib/gcc;/opt/slurm/slurm/lib64;/public4/soft/gcc/9.1.0/lib64;/lib64;/usr/lib64;/public4/soft/gcc/9.1.0/lib;/lib;/usr/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")




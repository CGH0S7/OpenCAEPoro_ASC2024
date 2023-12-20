set(CMAKE_CXX_COMPILER "/public4/soft/gcc/9.1.0/bin/c++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_VERSION "9.1.0")
set(CMAKE_CXX_PLATFORM_ID "Linux")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCXX 1)
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()




set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++;m;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/public4/soft/gcc/9.1.0/lib/gcc/x86_64-pc-linux-gnu/9.1.0;/public4/soft/gcc/9.1.0/lib/gcc;/opt/slurm/slurm/lib64;/public4/soft/gcc/9.1.0/lib64;/lib64;/usr/lib64;/opt/slurm/slurm/lib;/public4/soft/oneAPI/2022.1/mpi/2021.5.0/libfabric/lib;/public4/soft/oneAPI/2022.1/mpi/2021.5.0/lib/release;/public4/soft/oneAPI/2022.1/mpi/2021.5.0/lib;/public4/soft/oneAPI/2022.1/ccl/latest/lib/cpu;/public4/soft/oneAPI/2022.1/dal/latest/lib/intel64;/public4/soft/oneAPI/2022.1/tbb/latest/lib/intel64/gcc4.8;/public4/soft/oneAPI/2022.1/mkl/latest/lib/intel64;/public4/soft/oneAPI/2022.1/compiler/latest/linux/compiler/lib/intel64;/public4/soft/oneAPI/2022.1/compiler/latest/linux/lib;/public4/soft/oneAPI/2022.1/ipp/latest/lib;/public4/soft/gcc/9.1.0/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")




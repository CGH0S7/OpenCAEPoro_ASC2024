# Install script for directory: /home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/main

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro"
         RPATH "/usr/local/lib:/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/petsc_install/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro" TYPE EXECUTABLE FILES "/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/build/main/testOpenCAEPoro")
  if(EXISTS "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro"
         OLD_RPATH "/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/petsc_install/lib:::::::::::::::"
         NEW_RPATH "/usr/local/lib:/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/petsc_install/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/testOpenCAEPoro")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/build/main/CMakeFiles/testOpenCAEPoro.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
endif()

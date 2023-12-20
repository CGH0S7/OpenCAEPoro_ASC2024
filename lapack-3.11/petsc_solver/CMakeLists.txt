cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default configuration types and flags
set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Configs" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")

set(build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "${build_type}" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")

set(PETSC_DIR "/home/gh0s7/project/OpenCAEPoro_ASC2024/petsc-3.19.3/") 
set(PETSC_ARCH "petsc_install") 

project(petsc_solver)

include_directories(${PROJECT_SOURCE_DIR}/src/)
include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(${PETSC_DIR}/include/)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include/) #zhaoli, no modify
include_directories(/public1/home/sch10084/ParPennSim-4-14-16/common/include/)
include_directories(/public1/home/sch10084/ParPennSim-4-14-16/common/src/)
include_directories(/public1/home/sch10084/ParPennSim-4-14-16/Utility/include/)
include_directories(/public1/home/sch10084/ParPennSim-4-14-16/Utility/src/)

file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_library(petsc_solver STATIC ${SRC_FILES})

# Check for header file
find_path(PETSC_INCLUDE_DIRS petscksp.h
   HINTS ${PETSC_DIR}/include $ENV{PETSC_DIR}/include ${PROJECT_SOURCE_DIR}/petsc/include
   DOC "Directory where the PETSC header is located")
mark_as_advanced(PETSC_INCLUDE_DIRS)

# Check for petsc library
find_library(PETSC_LIBRARIES petsc
    HINTS ${PETSC_DIR}/${PETSC_ARCH}/lib $ENV{PETSC_DIR}/lib ${PROJECT_SOURCE_DIR}/PETSC/lib
    DOC "The PETSC library")
mark_as_advanced(PETSC_LIBRARIES)

set(PETSC_LIBRARIES ${PETSC_LIBRARIES})
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS})

# add_library(petsc_solver STATIC IMPORTED GLOBAL)
set_property(
    TARGET petsc_solver 
    PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
  set_property(
    TARGET petsc_solver 
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})

message(STATUS ${PETSC_INCLUDE_DIRS})
message(STATUS ${PETSC_LIBRARIES})
message(STATUS ${PETSC_DIR})

target_link_libraries(petsc_solver ${PETSC_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(petsc_solver
    "PETSC could not be found. Check PETSC_DIR."
    PETSC_LIBRARIES PETSC_INCLUDE_DIRS)


###############################################################################
## Installtion targets for lib and executable files
###############################################################################

# Set install location for petsc_solver library
set(CMAKE_INSTALL_LIBDIR "${PROJECT_SOURCE_DIR}/lib/")
install(TARGETS petsc_solver
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

## CMakeUninstall
#
# By default, CMake does not provide the "make uninstall" target. This is
# intended to prevent "make uninstall" from removing useful files from the
# system.
#
# This file provides the basis for a "uninstall" build target. It is used to
# delete the files listed in install_manifest.txt.
#
# This file is taken from
#   https://cmake.org/Wiki/CMake_FAQ#Can_I_do_.22make_uninstall.22_with_CMake.3F
#

if(NOT EXISTS "/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/build/install_manifest.txt")
  message(FATAL_ERROR "Cannot find install manifest: /home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/build/install_manifest.txt")
endif()

file(READ "/home/gh0s7/project/OpenCAEPoro_ASC2024/OpenCAEPoro/build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
foreach(file ${files})
  message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
  if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    exec_program(
      "/usr/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    if(NOT "${rm_retval}" STREQUAL 0)
      message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    endif()
  else(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  endif()
endforeach()

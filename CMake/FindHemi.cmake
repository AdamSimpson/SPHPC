# - Try to find IlmBase headers and library
# Once done this will define
#  HEMI_FOUND - System has hemi
#  HEMI_INCLUDE_DIRS - The hemi include directories
#  HEMI_DEFINITIONS - Compiler switches required for using hemi

find_package(PkgConfig)
pkg_check_modules(PC_HEMI Hemi)
set(HEMI_DEFINITIONS ${PC_HEMI_CFLAGS_OTHER})

find_path(HEMI_INCLUDE_DIR hemi/hemi.h
          HINTS ${PC_HEMI_INCLUDEDIR} ${PC_HEMI_INCLUDE_DIRS})

set(HEMI_INCLUDE_DIRS ${HEMI_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set HEMI to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Hemi  DEFAULT_MSG
                                  HEMI_INCLUDE_DIR)

mark_as_advanced(HEMI_INCLUDE_DIR)

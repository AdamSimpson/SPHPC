# - Try to find IlmBase headers and library
# Once done this will define
#  ILMBASE_FOUND - System has ilmbase
#  ILMBASE_INCLUDE_DIRS - The ilmbase include directories
#  ILMBASE_LIBRARIES - The libraries needed to use ilmbase
#  ILMBASE_DEFINITIONS - Compiler switches required for using ilmbase

find_package(PkgConfig)
pkg_check_modules(PC_ILMBASE IlmBase)
set(ILMBASE_DEFINITIONS ${PC_ILMBASE_CFLAGS_OTHER})

find_path(ILMBASE_INCLUDE_DIR OpenEXR/half.h
          HINTS ${PC_ILMBASE_INCLUDEDIR} ${PC_ILMBASE_INCLUDE_DIRS}
          PATH_SUFFIXES OpenEXR)

find_library(ILMBASE_LIBRARY NAMES Half Iex IlmThread Imath
             HINTS ${PC_ILMBASE_LIBDIR} ${PC_ILMBASE_LIBRARY_DIRS} )

set(ILMBASE_LIBRARIES ${ILMBASE_LIBRARY})
set(ILMBASE_INCLUDE_DIRS ${ILMBASE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ILMBASE to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(IlmBase  DEFAULT_MSG
                                  ILMBASE_LIBRARY ILMBASE_INCLUDE_DIR)

mark_as_advanced(ILMBASE_INCLUDE_DIR ILMBASE_LIBRARY)

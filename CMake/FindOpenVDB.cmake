# - Try to find OpenVDB headers and library
# Once done this will define
#  OPENVDB_FOUND - System has openvdb
#  OPENVDB_INCLUDE_DIRS - The openvdb include directories
#  OPENVDB_LIBRARIES - The libraries needed to use openvdb
#  OPENVDB_DEFINITIONS - Compiler switches required for using openvdb

# OpenVDB doesn't currently come with pkg-config functionality but perhaps someday
find_package(PkgConfig)
pkg_check_modules(PC_OPENVDB QUIET openvdb)
set(OPENVDB_DEFINITIONS ${PC_OPENVDB_CFLAGS_OTHER})

find_package(Ilmbase)
find_package(TBB)

find_path(OPENVDB_INCLUDE_DIR openvdb/openvdb.h
          HINTS ${PC_OPENVDB_INCLUDEDIR} ${PC_OPENVDB_INCLUDE_DIRS}
          PATH_SUFFIXES openvdb )

find_library(OPENVDB_LIBRARY NAMES openvdb
             HINTS ${PC_OPENVDB_LIBDIR} ${PC_OPENVDB_LIBRARY_DIRS} )

set(OPENVDB_LIBRARIES ${OPENVDB_LIBRARY} ${TBB_LIBRARIES} ${ILMBASE_LIBRARIES})
set(OPENVDB_INCLUDE_DIRS ${OPENVDB_INCLUDE_DIR} ${TBB_INCLUDE_DIR} ${ILMBASE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OPENVDB_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OpenVDB  DEFAULT_MSG
                                  OPENVDB_LIBRARY OPENVDB_INCLUDE_DIR)

mark_as_advanced(OPENVDB_INCLUDE_DIR OPENVDB_LIBRARY )

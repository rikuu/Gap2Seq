find_package(PkgConfig)
pkg_check_modules(PC_HTSLIB QUIET libhts)
set(HTSLIB_DEFINITIONS ${PC_LIBHTS_CFLAGS_OTHER})

find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
          HINTS ${PC_HTSLIB_INCLUDEDIR} ${PC_HTSLIB_INCLUDE_DIRS}
          PATH_SUFFIXES libhts)

find_library(HTSLIB_LIBRARY NAMES hts libhts
             HINTS ${PC_HTSLIB_LIBDIR} ${PC_HTSLIB_LIBRARY_DIRS})

get_filename_component(HTSLIB_LIBRARY_DIR ${HTSLIB_LIBRARY} DIRECTORY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSLIB DEFAULT_MSG
                                  HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR)

mark_as_advanced(HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)

set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
set(HTSLIB_LIBRARY_DIRS ${HTSLIB_LIBRARY_DIR})
set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})

#
# Find the native CGNS includes and library
#
# CGNS_INCLUDE_DIR - where to find cgns.h, etc.
# CGNS_LIBRARIES   - List of fully qualified libraries to link against when using CGNS.
# CGNS_FOUND       - Do not attempt to use CGNS if "no" or undefined.
find_path(CGNS_INCLUDE_DIR NAMES cgnslib_f.h
  /usr/local/include
  /usr/include
  PATHS ${PROJECT_SOURCE_DIR}/external/cgns/include
  )

find_library(CGNS_LIBRARY NAMES libcgns.a
  /usr/local/lib
  /usr/lib
  PATHS ${PROJECT_SOURCE_DIR}/external/cgns/lib
  )
# if cannot find static library in system libraries
# try to use external
if(CGNS_LIBRARY)
	message("CANNOT find libcgns.a in system libs, trying to use external/cgns")
	find_path(CGNS_INCLUDE_DIR NAMES cgnslib_f.h
	  PATHS ${PROJECT_SOURCE_DIR}/external/cgns/cgnslib_2.5
	)
	
	find_library(CGNS_LIBRARY NAMES libcgns.a
	  PATHS ${PROJECT_SOURCE_DIR}/external/cgns/cgnslib_2.5/LINUX
	)
endif(CGNS_LIBRARY)

set(CGNS_FOUND "NO")
if(CGNS_INCLUDE_DIR)
  if(CGNS_LIBRARY)
    set( CGNS_LIBRARIES ${CGNS_LIBRARY} )
    set( CGNS_FOUND "YES" )
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set CGNS_FOUND to TRUE if
# all listed variables are TRUE
#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(CGNS DEFAULT_MSG CGNS_LIBRARY CGNS_INCLUDE_DIR)


mark_as_advanced(
  CGNS_INCLUDE_DIR
  CGNS_LIBRARY
)

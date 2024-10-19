# If already in cache, be silent
if (SCOTCH_INCLUDES AND SCOTCH_LIBRARIES)
  set (SCOTCH_FIND_QUIETLY TRUE)
endif()

find_path(SCOTCH_INCLUDES NAMES ptscotch.h scotch.h HINTS ${EXTERNAL_LIB_DIR}/scotch/install/include)
find_library(SCOTCH_LIBRARIES NAMES libptscotcherr.a libptscotchparmetisv3.a libscotcherr.a libptscotch.a libptscotcherrexit.a libscotch.a libscotcherrexit.a  HINTS /usr/lib/scotch ${EXTERNAL_LIB_DIR}/scotch/install/lib)

if(NOT SCOTCH_DIR)
  get_filename_component(SCOTCH_DIR ${SCOTCH_INCLUDES} PATH)
  STRING(REGEX REPLACE "/include$" "" SCOTCH_DIR "${SCOTCH_DIR}")
  message(WARNING "${SCOTCH_DIR}")
  set (ENV{SCOTCH_DIR} "${SCOTCH_DIR}")
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH
                                  FOUND_VAR SCOTCH_FOUND
				  REQUIRED_VARS SCOTCH_INCLUDES SCOTCH_LIBRARIES SCOTCH_DIR
                                  VERSION_VAR SCOTCH_VERSION)
mark_as_advanced(SCOTCH_INCLUDES SCOTCH_LIBRARIES SCOTCH_DIR)

# If already in cache, be silent
IF(PUGIXML_INCLUDES AND PUGIXML_LIBRARIES)
  set(PUGIXML_FIND_QUIETLY TRUE)  
ENDIF()

find_path(PUGIXML_INCLUDES NAMES pugixml.hpp HINTS ${EXTERNAL_LIB_DIR}/pugixml-1.14/build/include)
set(PUGIXML_DIR "${PUGIXML_INCLUDES}/..")

find_library(PUGIXML_LIBRARIES NAMES libpugixml.a HINTS ${EXTERNAL_LIB_DIR}/pugixml-1.14/build/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PUGIXML FOUND_VAR PUGIXML_FOUND
                                  REQUIRED_VARS PUGIXML_INCLUDES PUGIXML_LIBRARIES PUGIXML_DIR)

mark_as_advanced(PUGIXML_INCLUDES PUGIXML_LIBRARIES PUGIXML_DIR)

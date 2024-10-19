# If already in cache, be silent
if (PASTIX_INCLUDES AND PASTIX_LIBRARIES)
  set (PASTIX_FIND_QUIETLY TRUE)
endif()

find_path(PASTIX_DIR NAMES PASTIXConfig.cmake  HINTS /usr/include/pastix ${EXTERNAL_LIB_DIR}/pastix/build)
find_path(PASTIX_INCLUDES NAMES pastix.h spm.h HINTS /usr/include/pastix ${EXTERNAL_LIB_DIR}/pastix/build/include)
find_library(PASTIX_LIBRARIES NAMES libpastix.a libpastix_kernels.a libspm.a HINTS /usr/lib/pastix ${EXTERNAL_LIB_DIR}/pastix/build/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PASTIX
                                  FOUND_VAR PASTIX_FOUND
                                  REQUIRED_VARS PASTIX_INCLUDES PASTIX_LIBRARIES PASTIX_DIR
                                  VERSION_VAR PASTIX_VERSION_STRING)

mark_as_advanced(PASTIX_DIR PASTIX_INCLUDES PASTIX_LIBRARIES PASTIX_VERSION)

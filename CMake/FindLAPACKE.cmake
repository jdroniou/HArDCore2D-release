if (LAPACKE_INCLUDE_DIRS)
  set (LAPACKE_FIND_QUIETLY TRUE)
endif()

find_path(LAPACKE_INCLUDE_DIRS NAMES lapacke.h HINTS ${EXTERNAL_LIB_DIR}/OpenBLAS/build/include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE
                                  FOUND_VAR LAPACKE_FOUND
                                  REQUIRED_VARS LAPACKE_INCLUDE_DIRS) 

mark_as_advanced(LAPACKE_INCLUDE_DIRS)








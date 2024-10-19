# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.


#set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH TRUE)
if (COMBBLAS_INCLUDES AND COMBBLAS_LIBRARIES)
  set(COMBBLAS_FIND_QUIETLY TRUE)
endif()

find_path(COMBBLAS_INCLUDES NAMES CombBLAS.h
	HINTS ${EXTERNAL_LIB_DIR}/CombBLAS/_build
	PATH_SUFFIXES include include/CombBLAS)
find_library(COMBBLAS_LIBRARIES NAMES CombBLAS
             HINTS ${EXTERNAL_LIB_DIR}/CombBLAS/_build
             PATH_SUFFIXES lib)

if(COMBBLAS_LIBRARIES AND COMBBLAS_INCLUDES)
  set(COMBBLAS_ROOT "${COMBBLAS_INCLUDES}/../..")
  set(COMBBLAS_INCS "-DTPL_COMBBLAS_INCLUDE_DIRS=${COMBBLAS_ROOT}/include\;${COMBBLAS_ROOT}/../Applications/BipartiteMatchings")
endif()

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args(COMBBLAS REQUIRED_VARS COMBBLAS_LIBRARIES COMBBLAS_INCLUDES COMBBLAS_INCS)

mark_as_advanced(COMBBLAS_INCLUDES COMBBLAS_LIBRARIES COMBBLAS_INCS)

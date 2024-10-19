# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

if (PARMETIS_INCLUDES AND PARMETIS_LIBRARIES)
  set(PARMETIS_FIND_QUIETLY TRUE)
endif (PARMETIS_INCLUDES AND PARMETIS_LIBRARIES)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".dylib" ".so")	
  set(LD_SEARCH $ENV{DYLD_LIBRARY_PATH})
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
  set(LD_SEARCH $ENV{LD_LIBRARY_PATH})
endif()

find_path(PARMETIS_INCLUDES NAMES parmetis.h PATHS
  ${EXTERNAL_LIB_DIR}/parmetis
  ENV CPATH
  PATH_SUFFIXES
  build/include
)

find_path(LOCAL_METIS_INCLUDES NAMES metis.h PATHS
  ${EXTERNAL_LIB_DIR}/parmetis
  ENV CPATH
  PATH_SUFFIXES
  build/include
  Linux-x86_64/metis/include
  build/Linux-x86_64/metis/include
)


find_library(PARMETIS_LIBRARIES NAMES parmetis 
  PATHS 
  ${EXTERNAL_LIB_DIR}/parmetis
  ENV LIBRARY_PATH 
  ${LD_SEARCH}
  ${LIB_INSTALL_DIR} 
  PATH_SUFFIXES build/lib lib)

find_library(LOCAL_METIS_LIBRARIES NAMES metis
  PATHS
  ${EXTERNAL_LIB_DIR}/parmetis
  ENV LIBRARY_PATH
  ${LD_SEARCH}
  ${LIB_INSTALL_DIR}
  PATH_SUFFIXES build/lib lib build/Linux-x86_64/libmetis Linux-x86_64/libmetis /build/Linux-x86_64/metis/build/lib/)

if(PARMETIS_LIBRARIES AND LOCAL_METIS_LIBRARIES AND PARMETIS_INCLUDES AND LOCAL_METIS_INCLUDES)
  set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARIES} ${LOCAL_METIS_LIBRARIES})
  set(PARMETIS_ROOT "${PARMETIS_INCLUDES}/..") # peut etre enlever " "
  set(PARMETIS_INCLUDES ${PARMETIS_INCLUDES} ${LOCAL_METIS_INCLUDES})
  set(PARMETIS_INCS "-DTPL_PARMETIS_INCLUDE_DIRS=${PARMETIS_ROOT}/include\;${PARMETIS_ROOT}/Linux-x86_64/metis/include")
  set(PARMETIS_LIBS "-DTPL_PARMETIS_LIBRARIES=${PARMETIS_ROOT}/lib/libparmetis.a\;${PARMETIS_ROOT}/Linux-x86_64/libmetis/libmetis.a")
endif()

find_package_handle_standard_args(ParMETIS
                                  REQUIRED_VARS PARMETIS_INCLUDES PARMETIS_LIBRARIES PARMETIS_ROOT PARMETIS_INCS PARMETIS_LIBS
)

mark_as_advanced(PARMETIS_INCLUDES PARMETIS_LIBRARIES PARMETIS_ROOT PARMETIS_INCS PARMETIS_LIBS)

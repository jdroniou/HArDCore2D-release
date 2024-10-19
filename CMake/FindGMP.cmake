
#set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH TRUE)
if (GMP_INCLUDES AND GMP_LIBRARIES)
  set(GMP_FIND_QUIETLY TRUE)
endif()

#
#
#    Hack in some MESSAGE calls in Modules/FindPkgConfig.cmake. In particular, before the execute_process call(s) if you want to dump how pkg-config is being invoked.
#
#    Now, like me, you might be astonished that your MESSAGE's never happen! What I learned is that if cmake invokes pkg-config for a particular package while checking dependencies for some other package earlier in the timeline, then the version it finds is cached in CMakeCache.txt which has session persistence. This means that if you later come and update that package (because some other dependency check is failing) then CMakeCache.txt still holds the old value and so despite you pulling your hair out, cmake errors out because the version or whatever isn't right.
#
#    Simply put, delete CMakeCache.txt if you update any packages as part of trying to get through cmake errors.

find_package(PkgConfig)
IF(PKGCONFIG_FOUND)  # first look for .pc in $ENV
  if(EXISTS "${EXTERNAL_LIB_DIR}/gmp-6.3.0/build/lib/pkgconfig") 
    set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${EXTERNAL_LIB_DIR}/gmp-6.3.0/build/lib/pkgconfig")
  endif()	  
  pkg_search_module(GMP "gmp>=6.1.2") # better stop at first successful match to avoid getting a list of directories
  if(GMP_FOUND)
    list(GET GMP_LINK_LIBRARIES 0 GMP_LIBRARY)
    set(GMP_INCLUDES ${GMP_INCLUDEDIR})
  endif()
  if(GMP_STATIC_FOUND)
    list(GET GMP_STATIC_LINK_LIBRARIES 0 GMP_STATIC)
    set(GMP_INCLUDES ${GMP_INCLUDEDIR})
  endif()
ENDIF()

if(NOT GMP_FOUND) # Switch to manual search
   find_path(GMP_INCLUDES NAMES gmp.h
             HINTS ${EXTERNAL_LIB_DIR}/gmp-6.3.0
	     PATH_SUFFIXES include Include build/Include build/include)
   find_library(GMP_LIBRARIES NAMES gmp
                HINTS ${EXTERNAL_LIB_DIR}/gmp-6.3.0
	        PATH_SUFFIXES lib build/lib)
   get_filename_component(GMP_LIBRARY ${GMP_LIBRARIES} REALPATH)
   string(REGEX MATCH "gmp-[0-9]+.[0-9]+.[0-9]+" GMP_LIB_VERSION ${GMP_LIBRARY}) # Retrieve version to check it match minimum required
   string(REGEX MATCH "-NOTFOUND" TMP ${GMP_LIBRARY})
   if((MPFR_LIB_VERSION STREQUAL "") OR (NOT TMP STREQUAL "")) # Version not found in path name, look in the .h
     file(STRINGS ${GMP_INCLUDES}/gmp.h GMP_VERSION_MAJOR REGEX "\#define __GNU_MP_VERSION ")
     file(STRINGS ${GMP_INCLUDES}/gmp.h GMP_VERSION_MINOR REGEX "\#define __GNU_MP_VERSION_MINOR")
     file(STRINGS ${GMP_INCLUDES}/gmp.h GMP_VERSION_PATCH REGEX "\#define __GNU_MP_VERSION_PATCHLEVEL")
     set(GMP_VERSION "${GMP_VERSION_MAJOR}.${GMP_VERSION_MINOR}.${GMP_VERSION_PATCH}")
   else()
     STRING(REGEX REPLACE "[\-a-zA-Z\n\(\)\ ]" "" GMP_VERSION "${GMP_LIB_VERSION}")
   endif()
endif()


include (FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP REQUIRED_VARS GMP_LIBRARIES GMP_INCLUDES VERSION_VAR GMP_VERSION)

mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES)

if (GMP_FOUND)
 message(STATUS "GMP version: ${GMP_VERSION}")
 message(STATUS "GMP include: ${GMP_INCLUDES}")
 message(STATUS "GMP library: ${GMP_LIBRARIES}")
endif()

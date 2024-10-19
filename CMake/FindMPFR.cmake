# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.


#set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH TRUE)
if (MPFR_INCLUDES AND MPFR_LIBRARIES)
  set(MPFR_FIND_QUIETLY TRUE)
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
  if(EXISTS "${EXTERNAL_LIB_DIR}/mpfr-4.2.1/build/lib/pkgconfig")
    set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${EXTERNAL_LIB_DIR}/mpfr-4.2.1/build/lib/pkgconfig")
  endif()	
  pkg_search_module(MPFR "mpfr>=4.0.2") # better stop at first successful match to avoid getting a list of directories
  if (MPFR_FOUND)
    list(GET MPFR_LINK_LIBRARIES 0 MPFR_LIBRARY)
    set(MPFR_INCLUDES ${MPFR_INCLUDEDIR})
  endif()
  if (MPFR_STATIC_FOUND)
    list(GET MPFR_STATIC_LINK_LIBRARIES 0 MPFR_STATIC)
    set(MPFR_INCLUDES ${MPFR_INCLUDEDIR})
  endif()
ENDIF()

if(NOT MPFR_FOUND) # Switch to manual search
   find_path(MPFR_INCLUDES NAMES mpfr.h
             HINTS ${EXTERNAL_LIB_DIR}/mpfr-4.2.1/build
             PATH_SUFFIXES include Include)
   find_library(MPFR_LIBRARIES NAMES mpfr
                HINTS ${EXTERNAL_LIB_DIR}/mpfr-4.2.1/build
	        PATH_SUFFIXES lib build)
   get_filename_component(MPFR_LIBRARY ${MPFR_LIBRARIES} REALPATH)
   string(REGEX MATCH "mpfr-[0-9]+.[0-9]+.[0-9]+" MPFR_LIB_VERSION ${MPFR_LIBRARY}) # Retrieve version to check it match minimum required
   if(MPFR_LIB_VERSION STREQUAL "") # Version not found in path name, look in the .h
     file(STRINGS ${MPFR_INCLUDES}/mpfr.h MPFR_VERSION REGEX "MPFR_VERSION_STRING")
   else()
     STRING(REGEX REPLACE "[\-a-zA-Z\n\(\)\ ]" "" MPFR_VERSION "${MPFR_LIB_VERSION}")
   endif()
endif()

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args(MPFR REQUIRED_VARS MPFR_LIBRARIES MPFR_INCLUDES VERSION_VAR MPFR_VERSION)

mark_as_advanced(MPFR_INCLUDES MPFR_LIBRARIES)

if (MPFR_FOUND)
 message(STATUS "MPFR version: ${MPFR_VERSION}")
 message(STATUS "MPFR include: ${MPFR_INCLUDES}")
 message(STATUS "MPFR library: ${MPFR_LIBRARIES}")
endif()

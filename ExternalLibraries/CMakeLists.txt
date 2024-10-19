#------------------------------------------------------------------------------#
#                             Helpful compilation libraries
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.5)  # to be commented
FIND_PACKAGE(Git QUIET)
if(GIT_FOUND)
 list(APPEND GIT_EXECUTABLE "clone")
endif()
FIND_PACKAGE(Wget QUIET)
if(NOT WGET_FOUND)
  FIND_PACKAGE(CURL)
  if(NOT CURL_FOUND)
    message(STATUS "Either curl or wget should be installed to download 3rd parts")
  else()
    set(GET_ONLINE_CMD "curl")
    list(APPEND GET_ONLINE_CMD "-LO")
    set(GET_ONLINE_CMD_FOUND CURL_FOUND)
  endif()
else()
  set(GET_ONLINE_CMD_FOUND WGET_FOUND)
  set(GET_ONLINE_CMD "wget")
endif()
set(PKGCONFIG_HINTS_PATH "$ENV{PATH}")
if(EXISTS "${EXTERNAL_LIB_DIR}/pkg-config-0.29.2/build/bin")
list(APPEND PKGCONFIG_HINTS_PATH "${EXTERNAL_LIB_DIR}/pkg-config-0.29.2/build/bin")
set(PKG_CONFIG_EXECUTABLE "${EXTERNAL_LIB_DIR}/pkg-config-0.29.2/build/pkg-config" )
endif()
find_program(PKG_EXE pkg-config HINTS ${PKGCONFIG_HINTS_PATH} NO_CACHE)
if(PKG_EXE)
 FIND_PACKAGE(PkgConfig REQUIRED QUIET)
else()
 set(PKGCONFIG_FOUND False)
endif()

function(clean_cache INCLUDES_PATH LIBRARIES_PATH) # Normally not required but safer
  unset(${INCLUDES_PATH}  CACHE)
  unset(${LIBRARIES_PATH} CACHE)
endfunction()

function(version_check LOC_EXE MIN_REQUIRED_VERSION ISFOUND)
# This function checks for version of a given executable (LOC_EXE) and compares it with a requirement (MIN_REQUIRED_VERSION)
# if conditions are filled it sets the variable (ISFOUND) at true. Note that ISFOUND is expected to be false at entry
# if the executable doesn't provide version, it is expected to exit with ISFOUND unchanged

# Transform output of `$ exe --version` in "major.minor.patch" format then convert it into cmake list
      execute_process(COMMAND ${LOC_EXE} --version OUTPUT_VARIABLE LOC_VERSION)
      string(REPLACE "\n" ";" TMP_LIST ${LOC_VERSION})
      list(GET TMP_LIST 0 LOC_VERSION)
      STRING(REGEX REPLACE "[\-a-zA-Z\n\(\)\ ]" "" LOC_VERSION "${LOC_VERSION}")
      string(REPLACE "." ";" VERSION_LIST ${LOC_VERSION})
      list(LENGTH VERSION_LIST VERSION_LENGTH)
# Convert MIN_REQUIRED_VERSION into cmake list
      string(REPLACE "." ";" MIN_REQUIRED_VERSION ${MIN_REQUIRED_VERSION})
      list(LENGTH MIN_REQUIRED_VERSION MIN_VERSION_LENGTH)

# Check if LOC_VERSION satifies MIN_REQUIRED_VERSION; ugly but does its job while adding arithmetic in cmake's list currently induces dependencies
      if ((MIN_VERSION_LENGTH GREATER_EQUAL 0) AND (VERSION_LENGTH GREATER_EQUAL 0))
        list(GET VERSION_LIST 0 LOC_VERSION_MAJOR)
        list(GET MIN_REQUIRED_VERSION 0 MIN_VERSION_MAJOR)
        # simply like found 3.xx.xx while requires 2.xx.xx
        if(LOC_VERSION_MAJOR GREATER MIN_VERSION_MAJOR)
            set(${ISFOUND} TRUE PARENT_SCOPE)
        endif()

        # case found 2.x.zz while requires 2.y.zz with x > y
        if( (LOC_VERSION_MAJOR EQUAL MIN_VERSION_MAJOR) AND (VERSION_LENGTH GREATER_EQUAL 1) )
           set(MIN_VERSION_MINOR 0)
           list(GET VERSION_LIST 1 LOC_VERSION_MINOR)
           list(GET MIN_REQUIRED_VERSION 1 MIN_VERSION_MINOR)

           if(LOC_VERSION_MINOR GREATER MIN_VERSION_MINOR)
              set(${ISFOUND} TRUE PARENT_SCOPE)
           endif()

           # case found 2.3.x while requires 2.3.y with x>y
           if( (LOC_VERSION_MINOR EQUAL MIN_VERSION_MINOR) AND (VERSION_LENGTH GREATER_EQUAL 2) )
               set(MIN_VERSION_PATCH 0)
               list(GET VERSION_LIST 2 LOC_VERSION_PATCH)
               list(GET MIN_REQUIRED_VERSION 2 MIN_VERSION_PATCH)

	       if(LOC_VERSION_PATCH GREATER_EQUAL MIN_VERSION_PATCH)
                   set(${ISFOUND} TRUE PARENT_SCOPE)
               endif()
           endif()
        endif()
      endif()
endfunction()

function(set_mpi_var PACK)
   set(_MPI_PKG ${PACK})
   pkg_check_modules(MPI_LANG_PKG "${_MPI_PKG}")
   if("${MPI_LANG_PKG_FOUND}")
     set(MPI_LANG_COMPILE_OPTIONS  ${MPI_LANG_PKG_CFLAGS}        CACHE STRING "MPI CXX compilation options"       FORCE PARENT_SCOPE)
     set(MPI_LANG_INCLUDE_PATH     ${MPI_LANG_PKG_INCLUDE_DIRS}  CACHE STRING "MPI CXX include directories"       FORCE PARENT_SCOPE)
     set(MPI_LANG_LINK_FLAGS       ${MPI_LANG_PKG_LDFLAGS}       CACHE STRING "MPI CXX linker flags"              FORCE PARENT_SCOPE)
     set(MPI_LANG_LIB_NAMES        ${MPI_LANG_PKG_LIBRARIES}     CACHE STRING "MPI CXX libraries to link against" FORCE PARENT_SCOPE)
     foreach(_MPI_LIB IN LISTS MPI_LANG_LIB_NAMES)
      if(_MPI_LIB)
        get_filename_component(_MPI_PLAIN_LIB_NAME "${_MPI_LIB}" NAME_WE   CACHE)
        get_filename_component(_MPI_LIB_NAME       "${_MPI_LIB}" NAME      CACHE)
        get_filename_component(_MPI_LIB_DIR        "${_MPI_LIB}" DIRECTORY CACHE)
        find_library(MPI_${_MPI_PLAIN_LIB_NAME}_LIBRARY NAMES "${_MPI_LIB_NAME}" "lib${_MPI_LIB_NAME}" HINTS ${_MPI_LIB_DIR} DOC "Location of the ${_MPI_PLAIN_LIB_NAME} library for MPI")
        mark_as_advanced(MPI_${_MPI_PLAIN_LIB_NAME}_LIBRARY)
      endif()
     endforeach()
   endif()
endfunction()

function(set_env_var BIN_PATH LIB_PATH INCLUDE_PATH PKG_CONFIG_PATH)
#~~~~ PKG path ~~~~#
if (NOT (PKG_CONFIG_PATH STREQUAL ""))
 set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}")
endif()
#~~~~ BIN PATHS ~~~~#
if (NOT (BIN_PATH STREQUAL ""))
 set(PATHV "$ENV{PATH}")
 if (NOT (PATHV STREQUAL ""))
  set(ENV{PATH} "$ENV{PATH}:${BIN_PATH}")
 else()
  set(ENV{PATH} "${BIN_PATH}")
 endif()
endif()
#~~~~ LIB PATHS ~~~~#
if (NOT (LIB_PATH STREQUAL ""))

 set(PATHV "$ENV{LIBRARY_PATH}")
 if (NOT (PATHV STREQUAL "") AND (NOT PATHV STREQUAL ":"))
  set(ENV{LIBRARY_PATH} "$ENV{LIBRARY_PATH}:${LIB_PATH}")
 else()
  set(ENV{LIBRARY_PATH} "${LIB_PATH}")
 endif()

 if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(PATHV "$ENV{DYLD_LIBRARY_PATH}")
  if (NOT (PATHV STREQUAL "") AND (NOT PATHV STREQUAL ":"))
   set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${LIB_PATH}")
  else()
   set(ENV{DYLD_LIBRARY_PATH} "${LIB_PATH}")
  endif()
 else()
  set(PATHV "$ENV{LD_LIBRARY_PATH}")
  if (NOT (PATHV STREQUAL "") AND (NOT PATHV STREQUAL ":"))
   set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${LIB_PATH}")
  else()
   set(ENV{LD_LIBRARY_PATH} "${LIB_PATH}")
  endif()
 endif()
endif()
#~~~~ INC PATHS ~~~~#
if (NOT (INCLUDE_PATH STREQUAL ""))
 set(PATHV "$ENV{CPATH}")
 if (NOT (PATHV STREQUAL "") AND (NOT PATHV STREQUAL ":"))
  set(ENV{CPATH} "$ENV{CPATH}:${INCLUDE_PATH}")
 else()
  set(ENV{CPATH} "${INCLUDE_PATH}")
 endif()
endif()
endfunction()

function(get_package getCMD getCMD_FOUND packageDIR packageSRC ArchFile BUILDDIR INSTALLDIR)
 file(GLOB RESULT ${packageDIR})
 #list(LENGTH RESULT RES_LEN)
 #IF ((EXISTS ${packageDIR}) AND (RES_LEN LESS_EQUAL 1))
 #  execute_process(COMMAND ${CMAKE_COMMAND} -E rm -rf ${packageDIR} COMMAND_ERROR_IS_FATAL ANY)
 # file(GLOB RESULT ${packageDIR})
 #ENDIF()
 IF (NOT EXISTS ${packageDIR})
   set(localGetCMD ${getCMD})
   list(APPEND localGetCMD "${packageSRC}")

   IF (getCMD_FOUND)
   message(STATUS "Downloading ...")
   execute_process(COMMAND ${localGetCMD} WORKING_DIRECTORY ${EXTERNAL_LIB_DIR} COMMAND_ERROR_IS_FATAL ANY)
   message(STATUS "Downloading done")
   IF (NOT ${ArchFile} STREQUAL "")
     message(STATUS "Extraction of the archive ...")
     file(ARCHIVE_EXTRACT INPUT ${EXTERNAL_LIB_DIR}/${ArchFile} DESTINATION ${EXTERNAL_LIB_DIR})
     if (NOT ${BUILDDIR} STREQUAL "")
      execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILDDIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     endif()
     execute_process(COMMAND ${CMAKE_COMMAND} -E rm -rf ${EXTERNAL_LIB_DIR}/${ArchFile} COMMAND_ERROR_IS_FATAL ANY)
     message(STATUS "Extraction of the archive done")
   endif()

   ELSE() # git/curl/wget not found
     message(FATAL_ERROR "requires ${getCMD} to download ${packageDIR}")
   ENDIF()
 ENDIF()
endfunction()

function(setFlex) # WARNING, later pass an installation directory to install for gcc then reinstall with gcc for pastix
# Installing flex if not present on the system (can be use by gcc or pastix)
# First look for flex in $PATH, it should work for /usr/bin, module env, ... ; if not found then look for local build
   set(FIND_EXE "flex")
   find_program(FLEX_EXE ${FIND_EXE} HINTS "$ENV{PATH}" NO_CACHE)
   set(FLEX_FOUND false)
   if (FLEX_EXE)
      version_check(${FLEX_EXE} "1.0.0" FLEX_FOUND)
      message(STATUS "FLEX 1")
   endif() # stop looking for flex in $PATH

# Then look for flex in local build directory
   if (NOT FLEX_FOUND) # 1
    unset(FLEX_EXE)
    find_program(FLEX_EXE ${FIND_EXE} PATHS ${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CACHE)
    if (FLEX_EXE)
      version_check(${FLEX_EXE} "1.0.0" FLEX_FOUND)
      if(FLEX_FOUND) # 2
        list(PREPEND APPEND_CMAKE_PREFIX_PATH "${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin")
        set(FLEX_BIN_DIR "${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin")
        set_env_var("${FLEX_BIN_DIR}" "${FLEX_BIN_DIR}/../lib" "${FLEX_BIN_DIR}/../include" "")
      endif()
    endif() # end FLEX_EXE FOUND
   endif() # sur (NOT FLEX_FOUND) 1

   IF(NOT FLEX_FOUND) # 3
    message(STATUS "flex was not found but is mandatory, try to fetch with wget")

# Creating variables to control installation directories
    SET(DIR_FLEX ${EXTERNAL_LIB_DIR}/flex-2.6.4)
    SET(BUILD_FLEX ${DIR_FLEX}/build)
    get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_FLEX} "https://github.com/westes/flex/files/981163/flex-2.6.4.tar.gz" "flex-2.6.4.tar.gz" ${BUILD_FLEX} "")

# runing configuration and installation
    message(STATUS "Configuring flex  ...")
    execute_process(COMMAND ./configure --prefix=${BUILD_FLEX} WORKING_DIRECTORY ${DIR_FLEX} RESULT_VARIABLE BUILD_RESULT )
    message(STATUS "Compiling flex ...")
    execute_process(COMMAND make         WORKING_DIRECTORY ${DIR_FLEX} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND make install WORKING_DIRECTORY ${DIR_FLEX} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    set_env_var("${BUILD_FLEX}/bin" "${BUILD_FLEX}/lib" "${BUILD_FLEX}/include" "")
    message(STATUS "Compiling flex done")
  endif() # sur FLEX NOT FOUND # 3
endfunction()

#-- ----------------------------------------------------------------------------#
#                             GCC if claimed
#-------------------------------------------------------------------------------#
if(INSTALL_GCC)
 SET(DIR_GCC ${EXTERNAL_LIB_DIR}/gcc)
 SET(BUILD_GCC ${EXTERNAL_LIB_DIR}/build_gcc)
 if(NOT EXISTS "${BUILD_GCC}")
 message(STATUS "Preparing GCC's installation")
 execute_process(COMMAND uname -m WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE ARCH_NAME)
 string(STRIP "${ARCH_NAME}" ARCH_NAME)
 if ((${CMAKE_SYSTEM_NAME} MATCHES "Darwin") AND (${ARCH_NAME} MATCHES "arm64"))
   # error is *** Configuration aarch64-apple-darwin23.3.0 not supported   Patch should only be for arm64 and removed as soon as gcc could, see :
   # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96168  https://unix.stackexchange.com/questions/730778/how-to-build-gcc-13-on-m1-macbook-pro
   get_package("${GIT_EXECUTABLE}" GIT_FOUND ${DIR_GCC} "https://github.com/iains/gcc-13-branch.git" "" "" "")
   execute_process(COMMAND curl -LO https://raw.githubusercontent.com/Homebrew/formula-patches/e923a0cd6c0e60bb388e8a5b8cd1dcf9c3bf7758/gcc/gcc-xcode15-warnings.diff WORKING_DIRECTORY "${EXTERNAL_LIB_DIR}/gcc-13-branch" COMMAND_ERROR_IS_FATAL ANY)
   #execute_process(COMMAND git apply gcc-xcode15-warnings.diff WORKING_DIRECTORY "${EXTERNAL_LIB_DIR}/gcc-13-branch" COMMAND_ERROR_IS_FATAL ANY)
   execute_process(COMMAND mv gcc-13-branch gcc WORKING_DIRECTORY ${EXTERNAL_LIB_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 else()
   get_package("${GIT_EXECUTABLE}" GIT_FOUND ${DIR_GCC} "https://gcc.gnu.org/git/gcc.git" "" "" "")
 endif()
 setFlex() # check if flex is found on the system
 execute_process(COMMAND ./contrib/download_prerequisites WORKING_DIRECTORY ${DIR_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)

 set(CONFIGURE_CMD ${DIR_GCC}/configure)
 if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  execute_process(COMMAND xcrun --show-sdk-path WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE SDKPATH)
  string(STRIP "${SDKPATH}" SDKPATH)
  set(ENV{SDKROOT} "${SDKPATH}")
  set(VARLD "DYLD_LIBRARY_PATH")
  if(NOT DEFINED ENV{CC}  AND (${CMAKE_C_COMPILER_ID}   MATCHES "[Cc]lang"))
    set(ENV{CC} "clang")
  endif()
  if(NOT DEFINED ENV{CXX} AND (${CMAKE_CXX_COMPILER_ID} MATCHES "[Cc]lang"))
    set(ENV{CXX} "clang++")
  endif()
  if(${ARCH_NAME} MATCHES "arm64")
   execute_process(COMMAND ${CONFIGURE_CMD} --prefix=${BUILD_GCC} --disable-multilib --enable-languages=all --with-system-zlib --with-ld=/Library/Developer/CommandLineTools/usr/bin/ld-classic --with-sysroot=${SDKPATH} WORKING_DIRECTORY ${BUILD_GCC} COMMAND_ERROR_IS_FATAL ANY)
  else()
   execute_process(COMMAND ${CONFIGURE_CMD} --prefix=${BUILD_GCC} --disable-multilib --enable-languages=all --with-sysroot=${SDKPATH} WORKING_DIRECTORY ${BUILD_GCC} COMMAND_ERROR_IS_FATAL ANY)
  endif()
  execute_process(COMMAND make -j 8 BOOT_CFLAGS="-O2 -Wno-error=deprecated-declarations" WORKING_DIRECTORY ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 else()
  execute_process(COMMAND ${CONFIGURE_CMD} --prefix=${BUILD_GCC} --disable-multilib --enable-languages=all WORKING_DIRECTORY ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
  execute_process(COMMAND make -j 8        WORKING_DIRECTORY ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
  set(VARLD "LD_LIBRARY_PATH")
 endif()	 
 execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
#----------- Configure and install -----------#
 message(STATUS "LIB: $ENV{LIBRARY_PATH}")
 #execute_process(COMMAND make -j 8        WORKING_DIRECTORY ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 execute_process(COMMAND make install     WORKING_DIRECTORY ${BUILD_GCC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 file(WRITE ${EXTERNAL_LIB_DIR}/env/set_gcc_env "export PATH=${BUILD_GCC}/bin:$PATH\nexport ${VARLD}=${BUILD_GCC}/lib64:$${BUILD_GCC}/lib:$${BUILD_GCC}:$${VARLD}\nexport CPATH=${BUILD_GCC}/include:$CPATH\n")
 endif() # BUILD_GCC Exists
#----------- Set environment variabales -----------#
 if(EXISTS "${BUILD_GCC}/bin/gcc-ranlib")
  set(LOC_RANLIB_EXE "${BUILD_GCC}/bin/gcc-ranlib")
 else()
  set(LOC_RANLIB_EXE "${BUILD_GCC}/bin/gcc-ar")
 endif() 
 list(APPEND CMAKE_PREFIX_PATH "${BUILD_GCC}")
 set_env_var("${BUILD_GCC}/bin" "${BUILD_GCC}/lib" "${BUILD_GCC}/include" "")
 set(ENV{CC}  "${BUILD_GCC}/bin/gcc")
 set(ENV{CXX} "${BUILD_GCC}/bin/g++")
 set(CMAKE_C_COMPILER          "${BUILD_GCC}/bin/gcc"     )
 set(CMAKE_CXX_COMPILER        "${BUILD_GCC}/bin/g++"     )
 set(CMAKE_Fortran_COMPILER    "${BUILD_GCC}/bin/gfortran")
 set(CMAKE_C_COMPILER_AR       "${BUILD_GCC}/bin/gcc-ar"  )
 set(CMAKE_C_COMPILER_RANLIB   "${LOC_RANLIB_EXE}"        )
 set(CMAKE_CXX_COMPILER_AR     "${BUILD_GCC}/bin/gcc-ar"  )
 set(CMAKE_CXX_COMPILER_RANLIB "${LOC_RANLIB_EXE}"        )
 unset(CMAKE_OSX_SYSROOT)
 set(CMAKE_C_COMPILER_ID       "GNU"                      )
 set(CMAKE_CXX_COMPILER_ID     "GNU"                      )
 set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_STD} ${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DOPTIMIZE -pthread -g" )
 set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -DOPTIMIZE -pthread -g" )
 set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -DOPTIMIZE -pthread -g" )
endif()	

#-- ----------------------------------------------------------------------------#
#                             Required libraries
#-------------------------------------------------------------------------------#
SET(Boost_NO_WARN_NEW_VERSIONS 1)
SET(Boost_DEBUG           "ON" CACHE_INTERNAL)
SET(Boost_USE_STATIC_LIBS "ON" CACHE_INTERNAL)

if(EXISTS "${EXTERNAL_LIB_DIR}/boost_1_84_0")
  list(APPEND CMAKE_PREFIX_PATH "${EXTERNAL_LIB_DIR}/boost_1_84_0")
endif()
FIND_PACKAGE(Boost COMPONENTS filesystem program_options timer chrono QUIET)
IF(NOT Boost_FOUND)
   unset(Boost_FOUND CACHE)
   clean_cache(Boost_DIR Boost_LIBRARIES)
   message(STATUS "Boost required components was not found, try to fetch with wget")
   SET(DIR_BOOST ${EXTERNAL_LIB_DIR}/boost_1_84_0/)
   SET(BUILD_BOOST ${DIR_BOOST})
   get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_BOOST} "https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz" "boost_1_84_0.tar.gz" ${BUILD_BOOST} "")
   message(STATUS "Compiling boost ...")
   if(${CMAKE_C_COMPILER_ID}   MATCHES "[Cc]lang")
    set(BOOST_TOOLSET "clang") #"darwin"
   else()
    set(BOOST_TOOLSET "gcc") 
    execute_process(COMMAND which ${CMAKE_CXX_COMPILER} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE COMPILER_PATH)
    string(STRIP "${COMPILER_PATH}" COMPILER_PATH)
    file(WRITE "${DIR_BOOST}/tools/build/src/user-config.jam" "using gcc : : ${COMPILER_PATH} ;")
   endif()
   execute_process(COMMAND ${DIR_BOOST}/bootstrap.sh --prefix=${BUILD_BOOST} --with-toolset=${BOOST_TOOLSET} WORKING_DIRECTORY ${DIR_BOOST} RESULT_VARIABLE BUILD_RESULT)
   message(STATUS "boostrap done, running b2...")
   execute_process(COMMAND ${DIR_BOOST}/b2 --prefix=${BUILD_BOOST} install   WORKING_DIRECTORY ${DIR_BOOST} RESULT_VARIABLE BUILD_RESULT)
   message(STATUS "b2 result: BUILD_RESULT")
   set(BOOST_ROOT "${BUILD_BOOST}" CACHE STRING INTERNAL)
ENDIF() # BOOST_NOT_FOUND

#------------------------------------------------------------------------------#
#      Checking for OpenMPI and installation if required and not found
#------------------------------------------------------------------------------#
if (USE_MPI)
  if(EXISTS "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build")
   set_env_var("${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/lib" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/include" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/lib/pkgconfig")
   SET(CMAKE_CXX_COMPILER     "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpic++" )
   SET(CMAKE_C_COMPILER       "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpicc"  )
   SET(CMAKE_Fortran_COMPILER "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpifort")
  endif()
  FIND_PACKAGE(MPI)
  if(NOT MPI_FOUND)
   message(STATUS "MPI was not found but you required it for scotch/pastix, try to fetch with wget")
# Creating variables to control installation directories
   SET(DIR_MPI ${EXTERNAL_LIB_DIR}/openmpi-5.0.0)
   SET(BUILD_MPI ${DIR_MPI}/build)
   get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_MPI} "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.0.tar.gz" "openmpi-5.0.0.tar.gz" ${BUILD_MPI} "")
   message(STATUS "Configuring MPI  ...")
   execute_process(COMMAND ./configure --prefix=${BUILD_MPI} WORKING_DIRECTORY ${DIR_MPI} RESULT_VARIABLE BUILD_RESULT)
   message(STATUS "Compiling MPI ...")
   execute_process(COMMAND make all install                  WORKING_DIRECTORY ${DIR_MPI} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)

# Set environment variables & paths for internal configuration
   set_env_var("${BUILD_MPI}/bin" "${BUILD_MPI}/lib" "${BUILD_MPI}/include" "${BUILD_MPI}/lib/pkgconfig")
   pkg_search_module(PKG_HWLOC REQUIRED IMPORTED_TARGET hwloc)
   SET(CMAKE_CXX_COMPILER     "${BUILD_MPI}/bin/mpic++" )
   SET(CMAKE_C_COMPILER       "${BUILD_MPI}/bin/mpicc"  )
   SET(CMAKE_Fortran_COMPILER "${BUILD_MPI}/bin/mpifort")

# Set environment variables & paths for compilation (make install) and later run (ex: mpirun ... ddr-magnetostatic ...) 
# NOTE: it need to be source such as: $ source [YOUR/LOCAL/PATH]/ExternalLibraries/env/set_mpi_env
   execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${EXTERNAL_LIB_DIR}/env" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
     set(VARLD "DYLD_LIBRARY_PATH")
   else()
     set(VARLD "LD_LIBRARY_PATH")
   endif()
   file(WRITE ${EXTERNAL_LIB_DIR}/env/set_mpi_env "export PATH=${BUILD_MPI}/bin:$PATH\nexport export ${VARLD}=${BUILD_MPI}/lib:$${VARLD}\nexport CPATH=${BUILD_MPI}/include:$CPATH\n")
   unset(VARLD)
   unset(HWLOC_NEED) # Internal use only
  endif() # MPI NOT FOUND via cmake's find_package
ENDIF() # USE_MPI


IF(USE_BLAS)
#--------------------------------------------------------------------------------#
# Looking for BLAS, if not found then try to look for openblas and install it
#--------------------------------------------------------------------------------#
 if(EXISTS "${EXTERNAL_LIB_DIR}/OpenBLAS/build/lib")
  set_env_var("${EXTERNAL_LIB_DIR}/OpenBLAS/build/bin" "${EXTERNAL_LIB_DIR}/OpenBLAS/build/lib" "${EXTERNAL_LIB_DIR}/OpenBLAS/build/include" "")
  list(PREPEND CMAKE_PREFIX_PATH ${EXTERNAL_LIB_DIR}/OpenBLAS)
 endif()
 find_package(BLAS)
 find_package(LAPACKE)
 if(BLAS_FOUND AND LAPACKE_FOUND)
    message(STATUS "BLAS found")
 else()
    unset(BLAS_FOUND CACHE)
    clean_cache(BLAS_INCLUDES BLAS_LIBRARIES)
    SET(DIR_BLAS ${EXTERNAL_LIB_DIR}/OpenBLAS)
    SET(BUILD_BLAS ${DIR_BLAS}/build)
    get_package("${GIT_EXECUTABLE}" GIT_FOUND ${DIR_BLAS} "https://github.com/OpenMathLib/OpenBLAS.git" "" "" "")
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_BLAS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND make         PREFIX=${BUILD_BLAS} WORKING_DIRECTORY ${DIR_BLAS} COMMAND_ERROR_IS_FATAL ANY) #OUTPUT_QUIET)
    execute_process(COMMAND make install PREFIX=${BUILD_BLAS} WORKING_DIRECTORY ${DIR_BLAS} COMMAND_ERROR_IS_FATAL ANY) #OUTPUT_QUIET)
    list(APPEND CMAKE_PREFIX_PATH ${DIR_BLAS})
    set_env_var("" "${BUILD_BLAS}/lib" "${BUILD_BLAS}/include" "${BUILD_BLAS}/lib/pkgconfig")
 endif() # BLAS_NOT_FOUND 
ENDIF() # USE_BLAS

#------------------------------------------------------------------------------#
#        Checking for Pastix and installation if required
#------------------------------------------------------------------------------#
IF(USE_PASTIX)
if(EXISTS "${EXTERNAL_LIB_DIR}/pastix/build/lib/cmake/pastix/")
  list(APPEND CMAKE_PREFIX_PATH ${EXTERNAL_LIB_DIR}/pastix/build/lib/cmake/pastix/)
endif()
FIND_PACKAGE(PASTIX)
IF(NOT PASTIX_FOUND)  # User's defined -DUSE_PASTIX=ON
  unset(PASTIX_FOUND CACHE)

#------------------------------------------------------------------------------------------------------------------------------------------------------------#	
# Pastix requires several dependencies for this implementation and/or to be build:
#  --- PACKAGE         (REQUIRED AT "build"==compilation / "current implementation"== compilation and run)
#  --- python > 3      (build)                                                                    NOTHING will be done if missing, user/admin install
#  --- mkl             (current implementation)                                                   NOTHING will be done if missing, user/admin install
#  --- openmpi         (current implementation)                                                   SOMETHING will be try if missing
#  --- hwloc           (build)  (need hwloc.pc in pkgconfig path)                                 SOMETHING will be try if missing
#   |--- pkg_config    (not provided, use `module load` or `sudo apt-get/ brew install /...`)     NOTHING will be done if missing, user/admin install
#  --- ptscotch        (current implementation)                                                   SOMETHING will be try if missing
#   |---  openmpi
#   |---  bison >= 3.4 (build)  (warning no version check perform by ptscotch, can fail at make)  SOMETHING will be try if missing
#   |---  git                                                                                     NOTHING will be done if missing, user/admin install
#   |---  flex         (build)                                                                    SOMETHING will be try if missing (minimum version unknown, set to 1.0.0 without objetive reason)
#   
#
#  if compilation fails for scotch due to unrealistic versions of parser_ll.c parser_ly.h parser_yy.c, be sure to use python > 3, then check for PYTHONPATH if it persists failing
#  as stantard system python could be `python2` and the executable `python3` is used for invocation of python > 3
#
#
#  Other requirements can be added if a given library don't perform well on a given cluster
#------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Useful package & environment variable
#------------------------------------------------------------------------------#
 string(REPLACE ":" ";" APPEND_CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")
#------------------------------------------------------------------------------#
#        Checking for Pastix and installation if required
#        Subcase pkg-config, mandatory for pastix to find hwloc 
#------------------------------------------------------------------------------#
if(NOT PKGCONFIG_FOUND)
 SET(DIR_PKG ${EXTERNAL_LIB_DIR}/pkg-config-0.29.2)
 SET(BUILD_PKG ${DIR_PKG}/build)
 get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_PKG} "http://pkgconfig.freedesktop.org/releases/pkg-config-0.29.2.tar.gz" "pkg-config-0.29.2.tar.gz" ${BUILD_PKG} "")
 message(STATUS "Configuring pkg-config  ...")
 execute_process(COMMAND ../configure --prefix=${BUILD_PKG} --with-internal-glib WORKING_DIRECTORY ${BUILD_PKG} RESULT_VARIABLE BUILD_RESULT)
 message(STATUS "Compiling pkg-config ...")
 execute_process(COMMAND make         WORKING_DIRECTORY ${BUILD_PKG} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_PKG} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 set_env_var("${BUILD_PKG}/bin" "" "" "")
 set(pkgconfig_DIR ${BUILD_PKG}/bin)
 unset(PKGCONFIG_FOUND)
 find_package(pkgconfig REQUIRED)
endif()

#------------------------------------------------------------------------------#
#        Subcase HWLOC (mpi dependency but still claimed by pastix-nompi)
#------------------------------------------------------------------------------#
# Hwloc is mandatory for pastix even if mpi is not use; if hwloc is not present, pastix's last step of installation will fail linking
# pastix/common/isched_nohwloc.c as no reference to pastix_warning function can be found (normal behavior: no such definition in pastix's sources)
# this behavior can be change if later versions of pastix (>6.3.0) are patched
   set(FIND_EXE "hwloc-ls")
   find_program(HWLOC_EXE ${FIND_EXE} HINTS "$ENV{PATH}" NO_CACHE)
   set(HWLOC_FOUND false)
   if (HWLOC_EXE)
      version_check(${HWLOC_EXE} "2.0.0" HWLOC_FOUND)
      if(HWLOC_FOUND)
        message(STATUS "hwloc found : ${HWLOC_EXE}")
        string(REPLACE ${FIND_EXE} "" HWLOC_BIN_DIR ${HWLOC_EXE})
        set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${HWLOC_BIN_DIR}/../lib/pkgconfig")	
        pkg_search_module(PKG_HWLOC REQUIRED IMPORTED_TARGET hwloc)
        set(BUILD_HWLOC "${HWLOC_BIN_DIR}../")
        list(APPEND APPEND_CMAKE_PREFIX_PATH "${BUILD_HWLOC}")
      endif()

   endif() # stop looking for hwloc in $PATH

# Then look for HWLOC in local build directory
   if (NOT HWLOC_FOUND)
    unset(HWLOC_EXE)
    find_program(HWLOC_EXE ${FIND_EXE} PATHS ${EXTERNAL_LIB_DIR}/hwloc-2.9.3/build/bin NO_SYSTEM_ENVIRONMENT_PATH NO_CACHE)
    if (HWLOC_EXE)
      version_check(${HWLOC_EXE} "2.0.0" HWLOC_FOUND)
      if(HWLOC_FOUND)
        string(REPLACE ${FIND_EXE} "" HWLOC_BIN_DIR ${HWLOC_EXE})
        set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${HWLOC_BIN_DIR}/../lib/pkgconfig")
        pkg_search_module(PKG_HWLOC REQUIRED IMPORTED_TARGET hwloc)
        list(APPEND APPEND_CMAKE_PREFIX_PATH "${BUILD_HWLOC}../")
      endif()
    endif() # end HWLOC_EXE FOUND
   endif() # stop looking for hwloc in local build directory


#  FIND_PACKAGE(hwloc HINTS ${EXTERNAL_LIB_DIR}/hwloc-2.9.3/)
  IF(NOT HWLOC_FOUND)
    unset(HWLOC_FOUND CACHE)
    message(STATUS "HWLOC was not found but mandatory for scotch/pastix, try to fetch with wget")
# Creating variables to control installation directories
     SET(DIR_HWLOC ${EXTERNAL_LIB_DIR}/hwloc-2.9.3)
     SET(BUILD_HWLOC ${DIR_HWLOC}/build)
     get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_HWLOC} "https://download.open-mpi.org/release/hwloc/v2.9/hwloc-2.9.3.tar.gz" "hwloc-2.9.3.tar.gz" ${BUILD_HWLOC} "")

# runing configuration and installation
     message(STATUS "Configuring hwloc  ...")
     execute_process(COMMAND ./configure --prefix=${BUILD_HWLOC} --enable-static   # should be tried with --enable-static=yes --disable-libnuma CFLAGS="-static"
                     WORKING_DIRECTORY ${DIR_HWLOC} RESULT_VARIABLE BUILD_RESULT)
     message(STATUS "Compiling hwloc ...")
     execute_process(COMMAND make                               WORKING_DIRECTORY ${DIR_HWLOC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)   
     execute_process(COMMAND make prefix=${BUILD_HWLOC} install WORKING_DIRECTORY ${DIR_HWLOC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     message(STATUS "Compilation done")

# exporting HWLOC environment variable for later use in cmake invocation by pastix (should be kept internal to this cmakelist)
      set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${BUILD_HWLOC}/lib/pkgconfig/")
# this looks for hwloc.pc file and creates a new target (RD RD hwloc is the new target)
# IMPORTED_TARGET requires CMake >= 3.6.3
      pkg_search_module(PKG_HWLOC REQUIRED IMPORTED_TARGET hwloc)  # (RD RD maybe give some HINTS and remove the set env)
      list(APPEND APPEND_CMAKE_PREFIX_PATH "${BUILD_HWLOC}")
    endif() # sur HWLOC NOT FOUND locally


#------------------------------------------------------------------------------#
#        Checking for Pastix and installation if required
#                Subcase BISON
#------------------------------------------------------------------------------#
# First look for bison in $PATH, it should work for /usr/bin, module env, ... ; if not found with suitable version (>2) then look for local build
   set(FIND_EXE "bison")
   find_program(LOC_EXE ${FIND_EXE} HINTS "$ENV{PATH}" NO_CACHE)
   set(BISON_FOUND false)
   if (LOC_EXE)
      version_check(${LOC_EXE} "3.0.0" BISON_FOUND)
      if(BISON_FOUND)
        set(BISON_EXE ${LOC_EXE})
      endif()
   endif() # stop looking for hwloc in $PATH

# Then look for bison in local build directory
   if (NOT BISON_FOUND)
    unset(LOC_EXE)
    find_program(LOC_EXE ${FIND_EXE} PATHS ${EXTERNAL_LIB_DIR}/bison-3.8.2/build/bin NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CACHE)
    if (LOC_EXE)
      version_check(${LOC_EXE} "3.0.0" BISON_FOUND)
      if(BISON_FOUND)
        set(BISON_EXE ${LOC_EXE})
      endif()
    endif() # end BISON_EXE FOUND
   endif() # sur (NOT BISON_FOUND) pour commencer a regarder dans local directory

   IF(NOT BISON_FOUND)
     unset(BISON_FOUND CACHE)
     message(STATUS "Bison was not found but mandatory for scotch/pastix, try to fetch with wget")
# Creating variables to control installation directories
     SET(DIR_BISON ${EXTERNAL_LIB_DIR}/bison-3.8.2)
     SET(BUILD_BISON ${DIR_BISON}/build)
     get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_BISON} "https://ftp.gnu.org/gnu/bison/bison-3.8.2.tar.gz" "bison-3.8.2.tar.gz" ${BUILD_BISON} "")

# runing configuration and installation
     message(STATUS "Configuring bison  ...")
     execute_process(COMMAND ./configure --prefix=${BUILD_BISON} WORKING_DIRECTORY ${DIR_BISON} RESULT_VARIABLE BUILD_RESULT)
     message(STATUS "Compiling bison ...")
     execute_process(COMMAND make         WORKING_DIRECTORY ${DIR_BISON} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make install WORKING_DIRECTORY ${DIR_BISON} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     set(BISON_EXE "${BUILD_BISON}/bin/bison")
     message(STATUS "Compiling done, bison is ${BISON_EXE}")
   endif() # sur BISON NOT FOUND

#------------------------------------------------------------------------------#
#        Checking for Pastix and installation if required
#                Subcase FLEX
#------------------------------------------------------------------------------#  
# Installing flex if not present on the system
# First look for flex in $PATH, it should work for /usr/bin, module env, ... ; if not found then look for local build
   set(FIND_EXE "flex")
   find_program(FLEX_EXE ${FIND_EXE} HINTS "$ENV{PATH}" NO_CACHE)
   set(FLEX_FOUND false)
   if (FLEX_EXE)
      version_check(${FLEX_EXE} "1.0.0" FLEX_FOUND)
   endif() # stop looking for flex in $PATH

# Then look for flex in local build directory
   if (NOT FLEX_FOUND)
    unset(FLEX_EXE)  
    find_program(FLEX_EXE ${FIND_EXE} PATHS ${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CACHE)
    if (FLEX_EXE)
      version_check(${FLEX_EXE} "1.0.0" FLEX_FOUND)
      if(FLEX_FOUND)
	list(PREPEND APPEND_CMAKE_PREFIX_PATH "${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin")
        set(FLEX_BIN_DIR "${EXTERNAL_LIB_DIR}/flex-2.6.4/build/bin") 
	set_env_var("${FLEX_BIN_DIR}" "${FLEX_BIN_DIR}/../lib" "${FLEX_BIN_DIR}/../include" "")
      endif()
    endif() # end FLEX_EXE FOUND
   endif() # sur (NOT FLEX_FOUND) pour commencer a regarder dans local directory

   IF(NOT FLEX_FOUND)
    unset(FLEX_FOUND CACHE)
    message(STATUS "flex was not found but mandatory for scotch/pastix, try to fetch with wget")

# Creating variables to control installation directories
    SET(DIR_FLEX ${EXTERNAL_LIB_DIR}/flex-2.6.4)
    SET(BUILD_FLEX ${DIR_FLEX}/build)
    get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_FLEX} "https://github.com/westes/flex/files/981163/flex-2.6.4.tar.gz" "flex-2.6.4.tar.gz" ${BUILD_FLEX} "")

# runing configuration and installation
    message(STATUS "Configuring flex  ...")
    execute_process(COMMAND ./configure --prefix=${BUILD_FLEX} WORKING_DIRECTORY ${DIR_FLEX} RESULT_VARIABLE BUILD_RESULT )
    message(STATUS "Compiling flex ...")
    execute_process(COMMAND make         WORKING_DIRECTORY ${DIR_FLEX} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND make install WORKING_DIRECTORY ${DIR_FLEX} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    set_env_var("${BUILD_FLEX}/bin" "${BUILD_FLEX}/lib" "${BUILD_FLEX}/include" "")
    message(STATUS "Compiling flex done")
  endif() # sur FLEX NOT FOUND

#------------------------------------------------------------------------------#
#        Checking for Pastix and installation if required
#                Subcase Scotch
#------------------------------------------------------------------------------#
#if(EXISTS "${EXTERNAL_LIB_DIR}/scotch/install/lib/cmake/scotch")
#  list(APPEND CMAKE_PREFIX_PATH ${EXTERNAL_LIB_DIR}/scotch/install/lib/cmake/scotch)
# endif()
 #FIND_PACKAGE(SCOTCH HINTS ${EXTERNAL_LIB_DIR}/scotch/install/lib/cmake/scotch)
 FIND_PACKAGE(SCOTCH)

 SET(INSTALL_SCOTCH ${DIR_SCOTCH}/install)
 IF(NOT SCOTCH_FOUND)
    unset(SCOTCH_FOUND CACHE)
    message(status "scotch not found")
    SET(DIR_SCOTCH ${EXTERNAL_LIB_DIR}/scotch)
    SET(BUILD_SCOTCH ${DIR_SCOTCH}/build)
    SET(INSTALL_SCOTCH ${DIR_SCOTCH}/install)
    get_package("${GIT_EXECUTABLE}" ${GIT_FOUND} ${DIR_SCOTCH} "https://gitlab.inria.fr/scotch/scotch.git" "" "" "")
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_SCOTCH}   OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${INSTALL_SCOTCH} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    set(ENV{CCS} CMAKE_CXX_COMPILER)
    set(ENV{CCD} CMAKE_CXX_COMPILER)
# configuration
      message(STATUS "Configuring scotch")
      if(USE_MPI)
        set(ENV{CCP} CMAKE_CXX_COMPILER)
        execute_process(COMMAND ${CMAKE_COMMAND} cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_SCOTCH} -DCMAKE_PREFIX_PATH=${APPEND_CMAKE_PREFIX_PATH} -DCMAKE_BUILD_TYPE=Release .. -DTHREADS=ON -DMPI_THREAD_MULTIPLE=ON -DBUILD_PTSCOTCH=ON -DCOMMON_PTHREAD_AFFINITY_LINUX=ON -DBISON_EXECUTABLE=${BISON_EXE}
                        WORKING_DIRECTORY ${BUILD_SCOTCH})
      else()
        execute_process(COMMAND ${CMAKE_COMMAND} cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_SCOTCH} -DCMAKE_PREFIX_PATH=${APPEND_CMAKE_PREFIX_PATH} -DCMAKE_BUILD_TYPE=Release .. -DTHREADS=ON -DMPI_THREAD_MULTIPLE=OFF -DBUILD_PTSCOTCH=OFF -DCOMMON_PTHREAD_AFFINITY_LINUX=ON -DBISON_EXECUTABLE=${BISON_EXE}
                        WORKING_DIRECTORY ${BUILD_SCOTCH})
      endif()
#-DCOMMON_OS_MACOS

# installation
      message(STATUS "Installing scotch")
      execute_process(COMMAND cmake --build ${BUILD_SCOTCH} --target install WORKING_DIRECTORY ${BUILD_SCOTCH})
      list(PREPEND APPEND_CMAKE_PREFIX_PATH "${INSTALL_SCOTCH}")
      set_env_var("" "${INSTALL_SCOTCH}/lib/" "${INSTALL_SCOTCH}/include/" "")
#   endif() # dir scotch not exists
ELSE() # SCOTCH IS FOUND
   list(PREPEND APPEND_CMAKE_PREFIX_PATH "${SCOTCH_INCLUDE_DIR}..")
ENDIF() #sur ptscotch not found


#------------------------------------------------------------------------------#
#       Pastix's installation
#------------------------------------------------------------------------------#
# dowload pastix with git
    SET(DIR_PASTIX ${EXTERNAL_LIB_DIR}/pastix)
    SET(BUILD_PASTIX ${DIR_PASTIX}/build)
    SET(PASTIX_SRC "-b;release-6.3.0;https://gitlab.inria.fr/solverstack/pastix.git")
    get_package("${GIT_EXECUTABLE}" GIT_FOUND ${DIR_PASTIX} "${PASTIX_SRC}" "" "" "")

# Setting directories variables and compilation environment variables
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_PASTIX} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    message(STATUS "Configuring pastix")
    if(USE_MPI)
      execute_process(COMMAND ${CMAKE_COMMAND} cmake .. -DCMAKE_INSTALL_PREFIX=${BUILD_PASTIX} -DCMAKE_PREFIX_PATH=${APPEND_CMAKE_PREFIX_PATH} -DPASTIX_INT64=OFF -DPASTIX_ORDERING_PTSCOTCH=ON -DPASTIX_WITH_MPI=ON -DPASTIX_WITH_FORTRAN=OFF
                      WORKING_DIRECTORY ${BUILD_PASTIX} COMMAND_ERROR_IS_FATAL ANY)
    else()
      execute_process(COMMAND ${CMAKE_COMMAND} cmake .. -DCMAKE_INSTALL_PREFIX=${BUILD_PASTIX} -DCMAKE_PREFIX_PATH=${APPEND_CMAKE_PREFIX_PATH} -DPASTIX_INT64=OFF -DPASTIX_ORDERING_SCOTCH=ON -DPASTIX_WITH_MPI=OFF -DPASTIX_WITH_FORTRAN=OFF
                     WORKING_DIRECTORY ${BUILD_PASTIX} COMMAND_ERROR_IS_FATAL ANY)
    endif()

# installation
   message(STATUS "Installing pastix")
   execute_process(COMMAND cmake --build . --target install WORKING_DIRECTORY ${BUILD_PASTIX})

# Pastix keeps fortran in its targets for cmake even when configured and built without   
   file(READ ${BUILD_PASTIX}/PASTIXConfig.cmake NEW_FILE_CONTENT)
   STRING(REPLACE "include(\"\$\{CMAKE_CURRENT_LIST_DIR\}/pastixfTargets.cmake\")" " " NEW_FILE_CONTENT "${NEW_FILE_CONTENT}" )  
   file(WRITE ${BUILD_PASTIX}/PASTIXConfig.cmake "${NEW_FILE_CONTENT}")
   file(READ ${BUILD_PASTIX}/lib/cmake/pastix/PASTIXConfig.cmake NEW_FILE_CONTENT)
   STRING(REPLACE "include(\"\$\{CMAKE_CURRENT_LIST_DIR\}/pastixfTargets.cmake\")" " " NEW_FILE_CONTENT "${NEW_FILE_CONTENT}" )
   file(WRITE ${BUILD_PASTIX}/lib/cmake/pastix/PASTIXConfig.cmake "${NEW_FILE_CONTENT}")
   clean_cache(PASTIX_INCLUDES PASTIX_LIBRARIES)
   list(APPEND CMAKE_PREFIX_PATH ${EXTERNAL_LIB_DIR}/pastix/build/lib/cmake/pastix/)
  endif() # PASTIX NOT FOUND
endif() # USE_PATIX

IF(USE_UMFPACK)
if(EXISTS "${EXTERNAL_LIB_DIR}/SuiteSparse/lib")
 list(PREPEND CMAKE_PREFIX_PATH ${EXTERNAL_LIB_DIR}/SuiteSparse)
endif()
FIND_PACKAGE(Umfpack)
IF(NOT UMFPACK_FOUND)
unset(UMFPACK_FOUND CACHE)
clean_cache(UMFPACK_INCLUDES UMFPACK_LIBRARIES)
#--------------------------------------------------------------------------------#
# UMFPACK has not been found, an attempt to built from SuiteSparse will be done
# It requires mpfr and gmp
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
# Looking for GMP > 6.1.2
#--------------------------------------------------------------------------------#
 find_package(GMP 6.1.2)
 if (NOT GMP_FOUND)
   unset(GMP_FOUND CACHE)
   SET(DIR_GMP ${EXTERNAL_LIB_DIR}/gmp-6.3.0)
   SET(BUILD_GMP ${DIR_GMP}/build)
   get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_GMP} "https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz" "gmp-6.3.0.tar.xz" ${BUILD_GMP} "")
   message(STATUS "Configuring GMP  ...")
   execute_process(COMMAND ./configure --prefix=${BUILD_GMP} WORKING_DIRECTORY ${DIR_GMP} RESULT_VARIABLE BUILD_RESULT )
   message(STATUS "Compiling GMP ...")
   execute_process(COMMAND make         WORKING_DIRECTORY ${DIR_GMP} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   execute_process(COMMAND make install WORKING_DIRECTORY ${DIR_GMP} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   set_env_var("" "${BUILD_GMP}/lib" "${BUILD_GMP}/include" "${BUILD_GMP}/lib/pkgconfig")
 endif()

#--------------------------------------------------------------------------------#
# Looking for MPFR > 4.0.2
#--------------------------------------------------------------------------------#
find_package(MPFR 4.0.2)
if (NOT MPFR_FOUND)
   unset(MPFR_FOUND CACHE)
   SET(DIR_MPFR ${EXTERNAL_LIB_DIR}/mpfr-4.2.1)
   SET(BUILD_MPFR ${DIR_MPFR}/build)
   get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_MPFR} "https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz" "mpfr-4.2.1.tar.xz" ${BUILD_MPFR} "")
   message(STATUS "Configuring MPFR  ...")
   execute_process(COMMAND ./configure --prefix=${BUILD_MPFR} WORKING_DIRECTORY ${DIR_MPFR} RESULT_VARIABLE BUILD_RESULT )
   message(STATUS "Compiling MPFR ...")
   execute_process(COMMAND make         WORKING_DIRECTORY ${DIR_MPFR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   execute_process(COMMAND make install WORKING_DIRECTORY ${DIR_MPFR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   set_env_var("" "${BUILD_MPFR}/lib" "${BUILD_MPFR}/include" "${BUILD_MPFR}/lib/pkgconfig")
 endif()

 message(STATUS "UMFPACK not found, try to fetch SuiteSparse with git")
 SET(DIR_SPARSE ${EXTERNAL_LIB_DIR}/SuiteSparse)
 #------------------------------------------------------------------------------#
 #  Checking if SuiteSparse has not been downloaded and extracted earlier
 #------------------------------------------------------------------------------#	

 get_package("${GIT_EXECUTABLE}" GIT_FOUND "${DIR_SPARSE}" "https://github.com/DrTimothyAldenDavis/SuiteSparse.git" "" "" "")

 message(STATUS "Building SuiteSparse  ...")
 set(BUILD_CMD "make")
 list(APPEND BUILD_CMD "local")
 execute_process(COMMAND ${BUILD_CMD} WORKING_DIRECTORY ${DIR_SPARSE} COMMAND_ERROR_IS_FATAL ANY)
 set(BUILD_CMD "make")
 list(APPEND BUILD_CMD "install")
 message(STATUS "Installing SuiteSparse  ...")
 execute_process(COMMAND ${BUILD_CMD} WORKING_DIRECTORY ${DIR_SPARSE} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
 set_env_var("${DIR_SPARSE}/bin" "${DIR_SPARSE}/lib" "${DIR_SPARSE}/include" "${DIR_SPARSE}/lib/pkgconfig")
 set(UMFPACK_INCS "${DIR_SPARSE}/include") 
 set(UMFPACK_LIBS "${DIR_SPARSE}/lib") 
# Set environment variables & paths for compilation (make install) and later run (ex: ./ddr-magnetostatic ...) 
# NOTE: it need to be source such as: $ source [YOUR/LOCAL/PATH]/ExternalLibraries/env/set_umfpack_env 
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${EXTERNAL_LIB_DIR}/env" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(VARLD "DYLD_LIBRARY_PATH")
  else()
    set(VARLD "LD_LIBRARY_PATH")
  endif()  
  file(WRITE ${EXTERNAL_LIB_DIR}/env/set_umfpack_env "export PATH=${DIR_SPARSE}/bin:$PATH\nexport ${VARLD}=${DIR_SPARSE}/lib:$${VARLD}\nexport CPATH=${DIR_SPARSE}/include:$CPATH\n")
endif() # UMFPACK NOT FOUND

ENDIF() # USE_UMFPACK
  
#------------------------------------------------------------------------------#
#  Checking if SuperLU is claimed then present and install it with it's
#  dependencies (BLAS, ) if needed.
#  Acitvation with -DUSE_SUPERLU=ON
#------------------------------------------------------------------------------#
IF(USE_SUPERLU)

if(EXISTS "${EXTERNAL_LIB_DIR}/superlu/bin")
 set_env_var("" "" "" "${EXTERNAL_LIB_DIR}/superlu/lib/pkgconfig")
endif()

if(EXISTS "${EXTERNAL_LIB_DIR}/superlu/lib")
  set_env_var("" "${EXTERNAL_LIB_DIR}/superlu/lib" "" "")
endif()
FIND_PACKAGE(SuperLU)
IF(NOT SUPERLU_FOUND)
unset(SUPERLU_FOUND CACHE)
clean_cache(SUPERLU_INCLUDES SUPERLU_LIBRARIES)
  if(USE_MPI) # NEED to check for Comblas and PARMETIS and install it if not found
  if(EXISTS "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build")
   set_env_var("${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/lib" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/include" "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/lib/pkgconfig")
   SET(CMAKE_CXX_COMPILER     "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpic++")
   SET(CMAKE_C_COMPILER       "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpicc")
   SET(CMAKE_Fortran_COMPILER "${EXTERNAL_LIB_DIR}/openmpi-5.0.0/build/bin/mpifort")
  endif()
  FIND_PACKAGE(MPI) # to be remove later and done separtely
#------------------------------------------------------------------------------#
#  Looking for CombBLAS
#------------------------------------------------------------------------------#
    find_package(COMBBLAS)
    if(NOT COMBBLAS_FOUND)
      unset(COMBBLAS_FOUND CACHE)
      message(STATUS "CombBLAS not found, try to fetch it with git")
      SET(DIR_COMBBLAS   ${EXTERNAL_LIB_DIR}/CombBLAS)
      SET(BUILD_COMBBLAS ${DIR_COMBBLAS}/_build)
      get_package("${GIT_EXECUTABLE}" GIT_FOUND ${DIR_COMBBLAS} "https://github.com/PASSIONLab/CombBLAS.git CombBLAS" "" "" "")  # WARNING
      execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_COMBBLAS}   OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
#------------------------------------------------------------------------------#
#  build CombBLAS
#------------------------------------------------------------------------------#
      execute_process(COMMAND ${CMAKE_COMMAND} cmake .. -DCMAKE_INSTALL_PREFIX=${BUILD_COMBBLAS}
                      WORKING_DIRECTORY ${BUILD_COMBBLAS} RESULT_VARIABLE BUILD_RESULT )
      execute_process(COMMAND make         WORKING_DIRECTORY ${BUILD_COMBBLAS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
      execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_COMBBLAS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
#------------------------------------------------------------------------------#
#  Set CombBLAS' variables to build superlu_dist
#------------------------------------------------------------------------------#
      set(COMBBLAS_ROOT ${BUILD_COMBBLAS})
      set(ENV{COMBBLAS_ROOT} ${BUILD_COMBBLAS})
      set(COMBBLAS_INCS "-DTPL_COMBBLAS_INCLUDE_DIRS=${COMBBLAS_ROOT}/include\;${COMBBLAS_ROOT}/../Applications/BipartiteMatchings")
      set_env_var("" "${BUILD_COMBBLAS}/lib" "${BUILD_COMBBLAS}/include" "")
      set_env_var("" "${BUILD_COMBBLAS}/lib" "${BUILD_COMBBLAS}/Applications/BipartiteMatchings" "")
    endif() # COMBLAS_NOT_FOUND
#------------------------------------------------------------------------------#
#  Looking for PARMETIS, MPI should be done earlier
#------------------------------------------------------------------------------#
    find_package(ParMETIS)
    if(NOT PARMETIS_FOUND)
      unset(PARMETIS_FOUND CACHE)
      message(STATUS "ParMETIS not found, try to fetch it with git")
      SET(DIR_PARMETIS   ${EXTERNAL_LIB_DIR}/parmetis)
      SET(BUILD_PARMETIS ${DIR_COMBBLAS}/build)
      set(parmetis_arch_file "parmetis-4.0.3.tar.gz")
      get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_COMBBLAS} "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/" ${parmetis_arch_file} "" "")
#------------------------------------------------------------------------------#
#  Extract and rename archive
#------------------------------------------------------------------------------#
#     file(ARCHIVE_EXTRACT INPUT ${EXTERNAL_LIB_DIR}/${parmetis_arch_file} DESTINATION ${EXTERNAL_LIB_DIR})
#     execute_process(COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/${parmetis_arch_file} COMMAND_ERROR_IS_FATAL ANY)
     file(GLOB TO_MOVE ${EXTERNAL_LIB_DIR}/parmetis*)
     execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${TO_MOVE} parmetis WORKING_DIRECTORY ${EXTERNAL_LIB_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_PARMETIS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/${parmetis_arch_file} COMMAND_ERROR_IS_FATAL ANY)
     message(STATUS "Extraction of the archive done")
#------------------------------------------------------------------------------#
#  build metis
#------------------------------------------------------------------------------#
     execute_process(COMMAND make config prefix=${DIR_PARMETIS}/metis/build WORKING_DIRECTORY ${DIR_PARMETIS}/metis OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make                                           WORKING_DIRECTORY ${DIR_PARMETIS}/metis OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make install                                   WORKING_DIRECTORY ${DIR_PARMETIS}/metis OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
#------------------------------------------------------------------------------#
#  build parmetis
#------------------------------------------------------------------------------#
     execute_process(COMMAND make config prefix=${DIR_PARMETIS}/build WORKING_DIRECTORY ${DIR_PARMETIS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make                                     WORKING_DIRECTORY ${DIR_PARMETIS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make install                             WORKING_DIRECTORY ${DIR_PARMETIS} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
#------------------------------------------------------------------------------#
#  Set metis & parmetis' variables to build superlu_dist
#------------------------------------------------------------------------------#
     set(ENV{PARMETIS_ROOT} ${DIR_PARMETIS}/build)
     set(PARMETIS_INCS "-DTPL_PARMETIS_INCLUDE_DIRS=$ENV{PARMETIS_ROOT}/include\;$ENV{PARMETIS_ROOT}/Linux-x86_64/metis/include")
     set(PARMETIS_LIBS "-DTPL_PARMETIS_LIBRARIES=$ENV{PARMETIS_ROOT}/lib/libparmetis.a\;$ENV{PARMETIS_ROOT}/Linux-x86_64/libmetis/libmetis.a")
     message(WARNING "PARMETIS_INCS ${PARMETIS_INCS}")
     message(WARNING "PARMETIS_LIBS ${PARMETIS_LIBS}")
     set(ENV{PARMETIS_BUILD} ${DIR_PARMETIS}/build/Linux-x86_64)
   endif() # PARMETIS_NOT_FOUND
  endif() # USE_MPI

 message(STATUS "SUPERLU not found, try to fetch it with git")
 SET(DIR_SUPERLU ${EXTERNAL_LIB_DIR}/superlu)
 SET(BUILD_SUPERLU ${DIR_SUPERLU}/build)
 # get_package(getCMD getCMD_FOUND packageDIR packageSRC ArchFile BUILDDIR INSTALLDIR)
 #------------------------------------------------------------------------------#
 #  Checking if SuperLU has not been downloaded and extracted earlier
 #------------------------------------------------------------------------------#
 file(GLOB RESULT ${DIR_SUPERLU})
 list(LENGTH RESULT RES_LEN)
 IF ((NOT EXISTS ${DIR_SUPERLU}) OR (RES_LEN EQUAL 0))
   IF (GIT_FOUND)
      if(USE_MPI)
        set(git_sources "https://github.com/xiaoyeli/superlu_dist.git")
      endif()
      if(SUPERLU_OPENMP AND (NOT USE_MPI))
	message(WARNING "SUPERLU_OPENMP")
        set(git_sources "https://github.com/xiaoyeli/superlu_mt.git")
      endif()
      if(SUPERLU_PTHREAD AND (NOT USE_MPI) AND (NOT SUPERLU_OPENMP))
        set(git_sources "https://github.com/xiaoyeli/superlu_mt.git")
      endif()
      if(NOT (USE_MPI OR SUPERLU_PTHREAD OR SUPERLU_OPENMP)) # use sequential
        set(git_sources "https://github.com/xiaoyeli/superlu.git")
      endif()
      execute_process(COMMAND "${GIT_EXECUTABLE}" ${git_sources} superlu WORKING_DIRECTORY ${EXTERNAL_LIB_DIR} COMMAND_ERROR_IS_FATAL ANY)
   ELSE() # git not found
     message(FATAL_ERROR "requires git to download SuperLU, BE SURE TO LOAD/INSTALL GIT, cmake failed successfully")
   ENDIF() # git found / not found
      if(USE_MPI)
        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
        set(FLAGS_C "-DCMAKE_C_FLAGS=-std=c99 -g -DPRNTlevel=0 -DDEBUGlevel=0")
        set(FLAGS_CXX "-std=c++11")
	set(MAKE_CMD "cmake")
	list(APPEND MAKE_CMD "${PARMETIS_INCS}")
	list(APPEND MAKE_CMD "${PARMETIS_LIBS}")
	list(APPEND MAKE_CMD "${COMBBLAS_INCS}")
	list(APPEND MAKE_CMD "${FLAGS_C}")

	execute_process(COMMAND ${MAKE_CMD} -DTPL_ENABLE_COMBBLASLIB=ON -DTPL_COMBBLAS_LIBRARIES=${COMBBLAS_ROOT}/lib/libCombBLAS.a -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=${FLAGS_CXX} -DTPL_ENABLE_BLASLIB=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${BUILD_SUPERLU} .. WORKING_DIRECTORY ${BUILD_SUPERLU} COMMAND_ECHO STDOUT)
        execute_process(COMMAND make         WORKING_DIRECTORY ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	# ensure eigen found libsuperlu when built
	execute_process(COMMAND ln -s libsuperlu_dist.a libsuperlu.a WORKING_DIRECTORY "${BUILD_SUPERLU}/lib" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
        set(SUPERLU_INCS "${BUILD_SUPERLU}/include") # to inform eigen of superlu's location when build
        set(SUPERLU_LIBS "${BUILD_SUPERLU}/lib")
      endif()
      if(SUPERLU_OPENMP AND (NOT USE_MPI))
        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${DIR_SUPERLU}/lib OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	file(COPY ${DIR_SUPERLU}/MAKE_INC/make.linux.openmp DESTINATION ${DIR_SUPERLU})
	file(RENAME ${DIR_SUPERLU}/make.linux.openmp ${DIR_SUPERLU}/make.inc)
	execute_process(COMMAND make blaslib WORKING_DIRECTORY ${DIR_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	execute_process(COMMAND make all     WORKING_DIRECTORY ${DIR_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	set_env_var("" "${DIR_SUPERLU}/lib" "" "")
	# ensure eigen found libsuperlu when built
	execute_process(COMMAND ln -s libsuperlu_mt_OPENMP.a libsuperlu.a WORKING_DIRECTORY "${DIR_SUPERLU}/lib" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	set(SUPERLU_INCS "${DIR_SUPERLU}/include") # to inform eigen of superlu's location when build
	set(SUPERLU_LIBS "${DIR_SUPERLU}/lib")
      endif()
      if(SUPERLU_PTHREAD AND (NOT USE_MPI) AND (NOT SUPERLU_OPENMP))
        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${DIR_SUPERLU}/lib OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	execute_process(COMMAND make blaslib WORKING_DIRECTORY ${DIR_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	execute_process(COMMAND make all     WORKING_DIRECTORY ${DIR_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	set_env_var("" "${DIR_SUPERLU}/lib" "" "")
	# ensure eigen found libsuperlu when built
	execute_process(COMMAND ln -s libsuperlu_mt_PTHREAD.a libsuperlu.a WORKING_DIRECTORY "${DIR_SUPERLU}/lib" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	set(SUPERLU_INCS "${DIR_SUPERLU}/include") # to inform eigen of superlu's location when build
	set(SUPERLU_LIBS "${DIR_SUPERLU}/lib")
      endif()
      if(NOT (USE_MPI OR SUPERLU_PTHREAD OR SUPERLU_OPENMP)) # use sequential
        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	execute_process(COMMAND ${CMAKE_COMMAND} cmake -DCMAKE_INSTALL_PREFIX=${BUILD_SUPERLU} ..
                        WORKING_DIRECTORY ${BUILD_SUPERLU})
        execute_process(COMMAND make         WORKING_DIRECTORY ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_SUPERLU} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
	set_env_var("" "" "" "${EXTERNAL_LIB_DIR}/superlu/lib/pkgconfig")
        set(SUPERLU_INCS "${BUILD_SUPERLU}/include") # to inform eigen of superlu's location when build
        set(SUPERLU_LIBS "${BUILD_SUPERLU}/lib")
      endif()
 ENDIF() # NOT EXISTS ${DIR_SUPERLU}
endif() # SUPERLU NOT FOUND

ENDIF()	# USE_SUPERLU

#------------------------------------------------------------------------------#
#        Looking for PugiXML if requested and installation if not found 
#------------------------------------------------------------------------------#
IF(USE_XML)
FIND_PACKAGE(PUGIXML)
  IF(NOT PUGIXML_FOUND)
     message(STATUS "PugiXML not found try to fetch with wget")
# Creating variables to control installation directories
     SET(DIR_PUGI ${EXTERNAL_LIB_DIR}/pugixml-1.14)
     SET(BUILD_PUGI ${DIR_PUGI}/build)
     get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_PUGI} "https://github.com/zeux/pugixml/releases/download/v1.14/pugixml-1.14.tar.gz" "pugixml-1.14.tar.gz" ${BUILD_PUGI} "")
# runing configuration and installation
     message(STATUS "Configuring PugiXML  ...")
     execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_PUGI} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND ${CMAKE_COMMAND} cmake -DCMAKE_INSTALL_PREFIX=${BUILD_PUGI} ..
                     WORKING_DIRECTORY ${BUILD_PUGI})
     message(STATUS "Compiling PugiXML ...")
     execute_process(COMMAND make         WORKING_DIRECTORY ${BUILD_PUGI} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_PUGI} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
     set_env_var("" "${BUILD_PUGI}/lib" "${BUILD_PUGI}/include" "${BUILD_PUGI}")
     message(STATUS "Compiling done")
     list(APPEND CMAKE_PREFIX_PATH ${BUILD_PUGI})
   ENDIF()
  unset(PUGIXML_INCLUDES CACHE)
  unset(PUGIXML_LIBRARIES CACHE)
  unset(PUGIXML_FOUND CACHE)
ENDIF()



#------------------------------------------------------------------------------#
#        Checking if Eigen is present. Eigen is a requirement of HArDCore
#------------------------------------------------------------------------------#
FIND_PACKAGE(Eigen3)
IF(NOT Eigen3_FOUND)
   SET(DIR_EIGEN ${EXTERNAL_LIB_DIR}/eigen-3.4.0)
   SET(BUILD_EIGEN ${DIR_EIGEN}/build)
   get_package("${GET_ONLINE_CMD}" ${GET_ONLINE_CMD_FOUND} ${DIR_EIGEN} "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" "eigen-3.4.0.tar.gz" ${BUILD_EIGEN} "")
   message(STATUS "Configuring Eigen  ...")
   set(EIGEN_MAKE_CMD "cmake")
   if(USE_UMFPACK)
     set(ENV{UMFPACKDIR} "${UMFPACK_INCS};${UMFPACK_LIBS}")
     set(ENV{CHOLMODDIR} "${UMFPACK_INCS};${UMFPACK_LIBS}")
   endif()
   if(USE_SUPERLU)
     set(ENV{SUPERLUDIR} "${SUPERLU_INCS};${SUPERLU_LIBS}")
     message(WARNING "SUPERLU eigen" $ENV{SUPERLUDIR})
   endif()

   if( (NOT ${CMAKE_C_COMPILER_ID} MATCHES "[Cc]lang") AND (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "[Cc]lang") )
     list(APPEND EIGEN_MAKE_CMD "-DCMAKE_C_FLAGS=-fopenmp -DCMAKE_CXX_FLAGS=-fopenmp")
     message(WARNING "ADDING -fopenmp")
   endif()
   execute_process(COMMAND ${EIGEN_MAKE_CMD} .. -DCMAKE_INSTALL_PREFIX=${BUILD_EIGEN} WORKING_DIRECTORY ${BUILD_EIGEN} RESULT_VARIABLE BUILD_RESULT )
   message(STATUS "Compiling Eigen ...")
   execute_process(COMMAND make install WORKING_DIRECTORY ${BUILD_EIGEN} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
   message(STATUS "Compiling done")
endif() # sur EIGEN3 NOT FOUND


################################################################################
#
# \file      cmake/FindMKL.cmake
# \author    Daniele A. Di Pietro
# \copyright 2020, Daniele A. Di Pietro, Université de Montpellier
# \brief     Find the Math Kernel Library from Intel
# \date      Thu 26 Jan 2017 02:05:50 PM MST
#
################################################################################

# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDES - MKL include files directories
#  MKL_RT - The MKL library
#
#  The environment variables MKLROOT and INTEL are used to find the library.
#  Everything else is ignored. If MKL is found "-DMKL_LP64" is added to
#  CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()

# If already in cache, be silent
if (MKL_INCLUDES AND MKL_LIBRARIES)
  set (MKL_FIND_QUIETLY TRUE)
endif()

if(NOT DEFINED ENV{MKLROOT})
  set (ENV{MKLROOT} "/usr")
endif()

find_path(MKL_INCLUDES 
  NAMES mkl.h 
  HINTS $ENV{MKLROOT}/include
        $ENV{MKLROOT}/include/mkl
        )

find_library(MKL_LIBRARIES
  NAMES mkl_rt
  PATHS $ENV{MKLROOT}/lib
  $ENV{MKLROOT}/lib/x86_64-linux-gnu
  $ENV{MKLROOT}/lib/intel64
  $ENV{INTEL}/mkl/lib/intel64
  NO_DEFAULT_PATH)

if (MKL_INCLUDE_DIR AND MKL_LIBRARIES)

  if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
      NOT DEFINED ENV{CRAY_PRGENVGNU} AND
      NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
      NOT DEFINED ENV{CRAY_PRGENVINTEL})
    set(ABI "-m64")
  endif()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_LP64 ${ABI}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_LP64 ${ABI}")
endif()

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDES)
MARK_AS_ADVANCED(MKL_INCLUDES MKL_LIBRARIES)

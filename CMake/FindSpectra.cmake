################################################################################
#
# \file      cmake/FindSpectra.cmake
# \author    Daniele A. Di Pietro, Jerome Droniou
# \copyright 2020, Daniele A. Di Pietro, Université de Montpellier; 2022, Jerome Droniou, Monash University
# \brief     Find the Spectra library
# \date      4 Dec 2022
#
################################################################################

# Find the Spectra Library
#
#  SPECTRA_FOUND - System has Spectra
#  SPECTRA_INCLUDES - Spectra include files directories
#
#  The environment variables SPECTRAROOT are used to find the library.
#  Everything else is ignored.
#
#  Example usage:
#
#  find_package(Spectra)
#  if(SPECTRA_FOUND)
#    target_link_libraries(TARGET ${SPECTRA_LIBRARIES})
#  endif()

# If already in cache, be silent
if (SPECTRA_INCLUDES AND SPECTRA_LIBRARIES)
  set (SPECTRA_FIND_QUIETLY TRUE)
endif()

find_path(SPECTRA_INCLUDES NAMES SymEigsBase.h HINTS /usr/include/Spectra include/Spectra)

# Handle the QUIETLY and REQUIRED arguments and set SPECTRA_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Spectra DEFAULT_MSG SPECTRA_INCLUDES)
MARK_AS_ADVANCED(SPECTRA_INCLUDES)



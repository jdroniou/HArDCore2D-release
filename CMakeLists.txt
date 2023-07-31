project(HArD::Core2D)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

#------------------------------------------------------------------------------#
#                             Compiler setup
#------------------------------------------------------------------------------#

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  SET(ADDITIONAL_FLAGS "-DWIN32")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#SET(CMAKE_CXX_COMPILER "c++")
SET(CMAKE_CXX_COMPILER "clang++")
SET(CMAKE_CXX_STD "-std=c++17")
#SET(CMAKE_CXX_FLAGS "-g -Wall -Wpedantic -Werror -Wno-uninitialized")
SET(CMAKE_CXX_FLAGS "-g -Wall -Wno-uninitialized")
SET(CMAKE_CXX_FLAGS_DEBUGALL "${CMAKE_CXX_STD} ${CMAKE_CXX_FLAGS} -DEIGEN_INITIALIZE_MATRICES_BY_NAN -fsanitize=address,undefined -pthread" )
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_STD} ${CMAKE_CXX_FLAGS} -O3 -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_STD} ${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_STD} ${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )

#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Debugall)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug Release RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
MESSAGE(STATUS "Build type : ${CMAKE_BUILD_TYPE}")

#------------------------------------------------------------------------------#
#                             Required libraries
#------------------------------------------------------------------------------#

SET(Boost_DEBUG "ON")
SET(Boost_USE_STATIC_LIBS "ON")
SET(HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(CMAKE_MODULE_PATH ${HOME_DIR}/CMake)

FIND_PACKAGE(Boost COMPONENTS filesystem program_options timer chrono REQUIRED)
SET(EXTERNAL_INCLUDES ${Boost_INCLUDE_DIR})
SET(EXTERNAL_LIBRARIES ${Boost_LIBRARIES})

FIND_PACKAGE(Eigen3 3.3.4 REQUIRED)
SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} "Eigen3::Eigen")
MESSAGE(STATUS "Eigen version: ${EIGEN3_VERSION_STRING} (${EIGEN3_INCLUDE_DIR})")

FIND_PACKAGE(SuperLU)
IF(SUPERLU_FOUND)
  SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${SUPERLU_LIBRARIES})
  SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${SUPERLU_INCLUDES})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_SUPERLU")
ENDIF()

FIND_PACKAGE(Umfpack)
IF(UMFPACK_FOUND)
  SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${UMFPACK_LIBRARIES})
  SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${UMFPACK_INCLUDES})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_UMFPACK")  
ENDIF()

FIND_PACKAGE(MKL)
IF(MKL_FOUND)
  SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${MKL_LIBRARIES})
  SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${MKL_INCLUDES})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MKL") 
ENDIF()

FIND_PACKAGE(Spectra)
IF(Spectra_FOUND)
  SET(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${SPECTRA_INCLUDE_DIR})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_SPECTRA")  
ENDIF()

INCLUDE_DIRECTORIES(${EXTERNAL_INCLUDES})

#------------------------------------------------------------------------------#
#                             Directories for compiled libraries
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(include/Math)
INCLUDE_DIRECTORIES(include/Mesh2D)
INCLUDE_DIRECTORIES(src/Mesh)
INCLUDE_DIRECTORIES(src/Quadrature)
INCLUDE_DIRECTORIES(src/Common)
INCLUDE_DIRECTORIES(src/HybridCore)
INCLUDE_DIRECTORIES(src/Plot)
INCLUDE_DIRECTORIES(src/DDRCore)
INCLUDE_DIRECTORIES(src/HHOCore)

ADD_SUBDIRECTORY(src/Mesh)
ADD_SUBDIRECTORY(src/Quadrature)
ADD_SUBDIRECTORY(src/Common)
ADD_SUBDIRECTORY(src/HybridCore)
ADD_SUBDIRECTORY(src/Plot)
ADD_SUBDIRECTORY(src/DDRCore)
ADD_SUBDIRECTORY(src/HHOCore)

#------------------------------------------------------------------------------#
#                             Directories for schemes
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(Schemes)
ADD_SUBDIRECTORY(Schemes)



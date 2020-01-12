project(HArD::Core2D)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
#------------------------------------------------------------------------------#
#                             Compiler setup
#------------------------------------------------------------------------------#

IF(CMAKE_SYSTEM_NAME=Windows)
  SET(ADDITIONAL_FLAGS "-DWIN32")
ENDIF(CMAKE_SYSTEM_NAME=Windows)

SET(CMAKE_CXX_COMPILER "c++")
SET(CMAKE_CXX_FLAGS_DEBUGALL "-std=c++14 -g -Wall -fsanitize=address,undefined" )
SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -Wall -O3 -DOPTIMIZE" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -Wall -O3 -DNDEBUG -DOPTIMIZE" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++14 -Wall -O3 -g -DNDEBUG -DOPTIMIZE" ${ADDITIONAL_FLAGS} )

#SET(CMAKE_BUILD_TYPE Debug)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build: Debug Release RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#------------------------------------------------------------------------------#
#                             Required libraries
#------------------------------------------------------------------------------#
SET(Boost_DEBUG "ON")
SET(Boost_USE_STATIC_LIBS "ON")
SET(HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(CMAKE_MODULE_PATH ${HOME_DIR}/CMake)

FIND_PACKAGE(Boost COMPONENTS filesystem program_options timer chrono REQUIRED)
FIND_PACKAGE(Eigen3 3.3.4)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

#------------------------------------------------------------------------------#
#                             Directories for compiled libraries
#------------------------------------------------------------------------------#
INCLUDE_DIRECTORIES(src/Mesh)
INCLUDE_DIRECTORIES(src/Quadrature)
INCLUDE_DIRECTORIES(src/HybridCore)
INCLUDE_DIRECTORIES(src/Plot)

ADD_SUBDIRECTORY(src/Mesh)
ADD_SUBDIRECTORY(src/Quadrature)
ADD_SUBDIRECTORY(src/HybridCore)
ADD_SUBDIRECTORY(src/Plot)

#------------------------------------------------------------------------------#
#                             Directories for schemes
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(Schemes)
ADD_SUBDIRECTORY(Schemes)


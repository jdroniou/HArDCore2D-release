project(HArD::Core2D)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

#------------------------------------------------------------------------------#
#                             Compiler setup
#------------------------------------------------------------------------------#

IF(CMAKE_SYSTEM_NAME=Windows)
  SET(ADDITIONAL_FLAGS "-DWIN32")
ENDIF(CMAKE_SYSTEM_NAME=Windows)

SET(CMAKE_CXX_COMPILER "c++")
SET(CMAKE_CXX_FLAGS "-std=c++17 -g -Wall -Wno-maybe-uninitialized")
SET(CMAKE_CXX_FLAGS_DEBUGALL "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined -pthread" )
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O3 -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DOPTIMIZE -pthread" ${ADDITIONAL_FLAGS} )

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
SET(EXTERNAL_INCLUDES ${Boost_INCLUDES})
SET(EXTERNAL_LIBRARIES ${Boost_LIBRARIES})

FIND_PACKAGE(Eigen3 3.3.4)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

#------------------------------------------------------------------------------#
#                             Directories for compiled libraries
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(src/Mesh)
INCLUDE_DIRECTORIES(src/Quadrature)
INCLUDE_DIRECTORIES(src/Common)
INCLUDE_DIRECTORIES(src/HybridCore)
INCLUDE_DIRECTORIES(src/Plot)

ADD_SUBDIRECTORY(src/Mesh)
ADD_SUBDIRECTORY(src/Quadrature)
ADD_SUBDIRECTORY(src/Common)
ADD_SUBDIRECTORY(src/HybridCore)
ADD_SUBDIRECTORY(src/Plot)

#------------------------------------------------------------------------------#
#                             Directories for schemes
#------------------------------------------------------------------------------#

INCLUDE_DIRECTORIES(Schemes)
ADD_SUBDIRECTORY(Schemes)


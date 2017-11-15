####################################################
#   Copyright (c) 2017 Idiap Research Institute    #
#   Written by James Newling <jnewling@idiap.ch>   #
####################################################

cmake_minimum_required(VERSION 3.0)

project(trimed)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(${PROJECT_NAME}_VERSION 0.1)
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 1)

option(BUILD_PYTHON_LIB "Build Python library" ON)
option(BUILD_SHARED_LIBS "Build shared library" ON)

set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" )

#Where to find FindCython.cmake
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

#find_package (Threads REQUIRED)

add_subdirectory(trimed)

if (BUILD_PYTHON_LIB)
  find_package(Cython REQUIRED)
  add_subdirectory(python)
endif()

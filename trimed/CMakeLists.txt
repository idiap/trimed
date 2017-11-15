#######################################################################
# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/   #
# Written by James Newling <jnewling@idiap.ch>                        #
#######################################################################

include(EnableCompilerWarnings)
add_definitions(-std=c++11)

include_directories("include")
file(GLOB_RECURSE source_files src/*cpp)
add_library("trimed" SHARED ${source_files})

# Inspiration from https://cmake.org/examples/
target_include_directories (trimed PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(trimed  ${CMAKE_THREAD_LIBS_INIT})

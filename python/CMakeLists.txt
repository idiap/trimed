#######################################################################
# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/   #
# Written by James Newling <jnewling@idiap.ch>                        #
#######################################################################

SET(CYTHONGEN_CPP ${CMAKE_CURRENT_BINARY_DIR}/cythongen/pytrimed.cpp)
SET(CYTHON_PYX ${CMAKE_CURRENT_SOURCE_DIR}/pytrimed.pyx)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cythongen)

include_directories("../trimed/include")

add_library(pytrimed SHARED ${CYTHONGEN_CPP})

SET_TARGET_PROPERTIES(pytrimed PROPERTIES LINKER_LANGUAGE CXX PREFIX "")

#How and when to cythonise
add_custom_command(
  OUTPUT ${CYTHONGEN_CPP}
  DEPENDS ${CYTHON_PYX}
  COMMAND cython --cplus ${CYTHON_PYX} -o ${CYTHONGEN_CPP}
)

target_include_directories(pytrimed PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(pytrimed LINK_PUBLIC trimed)

execute_process(COMMAND python-config --cflags OUTPUT_VARIABLE CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REPLACE "-Wall" "" CFLAGS ${CFLAGS})
string(REPLACE "-Wextra" "" CFLAGS ${CFLAGS})
string(REPLACE "-pedantic" "" CFLAGS ${CFLAGS})

SET(CMAKE_CXX_FLAGS  " ${CFLAGS} --std=c++11")


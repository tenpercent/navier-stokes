cmake_minimum_required (VERSION 2.6)

MESSAGE (STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR})
MESSAGE (STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR})

IF (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  MESSAGE (SEND_ERROR "you should build in ./build directory")
  MESSAGE (FATAL_ERROR "remove CMakeCache.txt manually to rebuild")
ENDIF ()

set (CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -pedantic -std=c99 -O2")
set (CMAKE_C_FLAGS_DEBUG "-g -pg -Wall -Wextra -pedantic -std=c99 -O0")

option (PDESolver_NO_LASPack "Build without LASPack" OFF)
option (PDESolver_COLORFUL_OUTPUT "Use fancy output style" OFF)

project (PDESolver)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory (src)

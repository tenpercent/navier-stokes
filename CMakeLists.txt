cmake_minimum_required (VERSION 2.6)

MESSAGE (STATUS "CMAKE_BINARY_DIR:         " ${CMAKE_BINARY_DIR})
MESSAGE (STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR})

IF (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  MESSAGE (SEND_ERROR "you should build in ./build directory")
  MESSAGE (FATAL_ERROR "remove CMakeCache.txt manually to rebuild")
ENDIF ()

project (PDESolver C)

option (PDESolver_NO_LASPack "Build without LASPack" OFF)
option (PDESolver_COLORFUL_OUTPUT "Use fancy output style" OFF)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic --std=c90")

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory (laspack-1.12.2)
add_subdirectory (pdesolver)

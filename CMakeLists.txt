cmake_minimum_required (VERSION 2.6)

set (CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -pedantic -std=c99 -O2")
set (CMAKE_C_FLAGS_DEBUG "-g -pg -Wall -Wextra -pedantic -std=c99 -O0")

project (PDESolver)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory (src)

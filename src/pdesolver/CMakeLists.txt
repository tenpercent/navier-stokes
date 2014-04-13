include_directories(${PROJECT_SOURCE_DIR}/include/pdesolver)
include_directories(${PROJECT_SOURCE_DIR}/include/laspack)
include_directories(${PROJECT_SOURCE_DIR}/include)

set (PDESolver_sources
  main.c
  initialize.c
  construction.c
  calculation.c
  norm.c
  functions.c
  iterative_method.c
  sparse_matrix.c
  print.c
  export.c
)

add_executable (PDESolver ${PDESolver_sources})

target_link_libraries (PDESolver LASPack)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
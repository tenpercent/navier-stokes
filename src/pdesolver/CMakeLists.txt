include_directories(${PROJECT_SOURCE_DIR}/include/pdesolver)

if (NOT PDESolver_NO_LASPack)
  include_directories(${PROJECT_SOURCE_DIR}/include/laspack)
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif (NOT PDESolver_NO_LASPack)

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

if (PDESolver_NO_LASPack)
  add_definitions(-DNO_LASPACK)
endif (PDESolver_NO_LASPack)

if (PDESolver_COLORFUL_OUTPUT)
  add_definitions(-DALTERNATIVE_OUTPUT)
endif (PDESolver_COLORFUL_OUTPUT)

add_executable (PDESolver ${PDESolver_sources})

if (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver LASPack m)
else (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver m)
endif (NOT PDESolver_NO_LASPack)

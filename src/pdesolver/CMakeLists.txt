include_directories(${PROJECT_SOURCE_DIR}/include/pdesolver)

if (NOT PDESolver_NO_LASPack)
  include_directories(${PROJECT_SOURCE_DIR}/include/laspack)
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif (NOT PDESolver_NO_LASPack)

set (PDESolver_smooth_sources
  pdesolver.c
  initialize.c
  construction.c
  calculation.c
  fill_system_smooth.c
  norm.c
  functions.c
  iterative_method.c
  sparse_matrix.c
  print.c
  export.c
)

set (PDESolver_non_smooth_sources
  pdesolver.c
  initialize.c
  construction.c
  calculation.c
  fill_system_non_smooth.c
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

add_executable (PDESolver_12 ${PDESolver_smooth_sources})
add_executable (PDESolver_34 ${PDESolver_non_smooth_sources})

if (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver_12 LASPack m)
  target_link_libraries (PDESolver_34 LASPack m)
else (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver_12 m)
  target_link_libraries (PDESolver_34 m)
endif (NOT PDESolver_NO_LASPack)



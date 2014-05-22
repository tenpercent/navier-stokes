include_directories(${PROJECT_SOURCE_DIR}/pdesolver/include/)

if (NOT PDESolver_NO_LASPack)
  include_directories(${PROJECT_SOURCE_DIR}/laspack-1.12.2)
endif (NOT PDESolver_NO_LASPack)

set (PDESolver_common_sources 
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

set (PDESolver_12_specific_sources
  start_conditions_12.c
  fill_system_smooth.c
  pdesolver_smooth.c
)

set (PDESolver_3_specific_sources
  start_conditions_3.c
  fill_system_non_smooth.c
  pdesolver_non_smooth.c  
)

set (PDESolver_4_specific_sources
  fill_system_non_smooth.c
  start_conditions_4.c
  pdesolver_non_smooth.c 
)

if (PDESolver_NO_LASPack)
  add_definitions(-DNO_LASPACK)
endif (PDESolver_NO_LASPack)

if (PDESolver_COLORFUL_OUTPUT)
  add_definitions(-DALTERNATIVE_OUTPUT)
endif (PDESolver_COLORFUL_OUTPUT)

add_executable (PDESolver_12 ${PDESolver_common_sources} 
                             ${PDESolver_12_specific_sources})
add_executable (PDESolver_3 ${PDESolver_common_sources} 
                            ${PDESolver_3_specific_sources}
)
add_executable (PDESolver_4 ${PDESolver_common_sources} 
                            ${PDESolver_4_specific_sources}
)

if (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver_12 laspack m)
  target_link_libraries (PDESolver_3 laspack m)
  target_link_libraries (PDESolver_4 laspack m)
else (NOT PDESolver_NO_LASPack)
  target_link_libraries (PDESolver_12 m)
  target_link_libraries (PDESolver_3 m)
  target_link_libraries (PDESolver_4 m)
endif (NOT PDESolver_NO_LASPack)



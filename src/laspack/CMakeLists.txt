include_directories(${PROJECT_SOURCE_DIR}/include)

set (LASPack_sources
  eigenval.c
  errhandl.c
  factor.c
  itersolv.c
  matrix.c
  mlsolv.c
  operats.c
  precond.c
  qmatrix.c
  rtc.c
  vector.c
)

add_library (LASPack ${LASPack_sources})
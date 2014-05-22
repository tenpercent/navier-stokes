include_directories(..)

set(laspack_sources
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
  vector.c)

add_library (laspack ${laspack_sources})

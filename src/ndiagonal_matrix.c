#include "ndiagonal_mathix.h"
#include <stdlib.h>

#define INDEX(i,j) ()

typedef struct {
  double * p_to_begin;
  unsigned side;
  unsigned max_diagonal;
  unsigned min_diagonal;

  (void *) (*allocate) (unsigned);
  void (*deallocate) (void *);
} Almost_diagonal_matrix;

void Almost_diagonal_matrix_Construct (Almost_diagonal_matrix * matrix,
    unsigned side,
    unsigned max_diagonal,
    unsigned min_diagonal) {

  matrix->allocate = malloc;
  matrix->deallocate = free;

  matrix->side = side;
  matrix->max_diagonal = max_diagonal;
  matrix->min_diagonal = min_diagonal;

  matrix->p_to_begin = (double *) matrix->allocate(
    (max_diagonal * side - 
    (max_diagonal - min_diagonal) * (max_diagonal - min_diagonal + 1)) 
      * sizeof (double));
  return;
}

void Almost_diagonal_matrix_Destruct (Almost_diagonal_matrix * matrix) {
  matrix->deallocate(matrix->p_to_begin);
  return;
}

void Almost_diagonal_matrix_Apply_to_vector (
    Almost_diagonal_matrix * matrix,
    double * source_vector,
    double * result_vector) {

  unsigned incomplete_rows = matrix->max_diagonal - matrix->min_diagonal;

  for (unsigned i = 0; i < incomplete_rows; ++i) {
    result_vector[i] = 0;
    for (unsigned j = 0; j < incomplete_rows + i; ++j) {
      result_vector[i] += source_vector[j] * matrix[]
    }
  }
}

#undef INDEX(i,j)

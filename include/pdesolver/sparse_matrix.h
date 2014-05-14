#ifndef _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H

typedef struct {
  /* public stuff */
  unsigned size;
  unsigned *indices;
  double *elements;
  /* for private use */
  unsigned current_row;
  unsigned filled_element;
} Sparse_matrix;

void Sparse_matrix_Construct (
    Sparse_matrix * this,
    unsigned size,
    unsigned nzcount);

void Sparse_matrix_Clear (
    Sparse_matrix * this);

void Sparse_matrix_Append_element (
    Sparse_matrix * this,
    unsigned row,
    unsigned col,
    double value);

void Sparse_matrix_Finish_filling (
    Sparse_matrix * this);

void Sparse_matrix_Destruct (
    Sparse_matrix * this);

void Sparse_matrix_Apply_to_vector (
    Sparse_matrix const * this,
    double const * vector,
    double * newVector);

#ifndef NO_LASPACK
#include <qmatrix.h>
void Sparse_matrix_to_QMatrix (
    Sparse_matrix const * this,
    QMatrix * qmatrix);
#endif /* NO_LASPACK */

#endif /* _SPARSE_MATRIX_H */

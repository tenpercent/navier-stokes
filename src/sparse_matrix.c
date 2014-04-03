#include <stdlib.h>
#include "sparse_matrix.h"

/* Useful for arrays construction */
#define NEW(type, len) ((type*)malloc((len)*sizeof(type)))

void Sparse_matrix_Construct (
    Sparse_matrix * this,
    unsigned size,
    unsigned nzcount) {

  this->size = size;
  this->indices = NEW(unsigned, size + nzcount + 1);
  this->elements = NEW(double, size + nzcount + 1);

  Sparse_matrix_Clear (this);
}

void Sparse_matrix_Clear (
    Sparse_matrix * this) {

  this->current_row = 0;
  this->indices[0] = this->size + 1;
  this->filled_element = this->size;
}

/* We store elements in Modified Compressed Sparse Row
 * (MSR) format */

void Sparse_matrix_Append_element (
    Sparse_matrix * this,
    unsigned row,
    unsigned col,
    double value) {

  if (row == col) {
    this->elements[row] = value;
    return;
  }

  ++(this->filled_element);
  while (this->current_row < row) {
    ++(this->current_row);
    this->indices[this->current_row] = this->filled_element;
  }
  this->indices[this->filled_element] = col;
  this->elements[this->filled_element] = value;
}

void Sparse_matrix_Finish_filling (
    Sparse_matrix * this) {

  ++this->filled_element;
  while (this->current_row < this->size) {
    ++(this->current_row);
    this->indices[this->current_row] = this->filled_element;
  }
}

void Sparse_matrix_Destruct (
    Sparse_matrix * this) {

  free(this->indices);
  free(this->elements);
}

void Sparse_matrix_Apply_to_vector (
    Sparse_matrix * this,
    double * vector,
    double * newVector) {

  unsigned i, ind;
  for (i = 0; i < this->size; ++i) {
    newVector[i] = this->elements[i] * vector[i];
    for (ind = this->indices[i]; ind < this->indices[i+1]; ++ind) {
      newVector[i] += this->elements[ind] * vector[this->indices[ind]];
    }
  }
}

#ifndef NO_LASPACK
void Sparse_matrix_to_QMatrix (
    Sparse_matrix const * this,
    QMatrix * qmatrix) {

  unsigned row, index, row_length, position;
  for (row = 0; row < this->size; ++row) {
    row_length = this->indices[row + 1] - this->indices[row];
    Q_SetLen (qmatrix, row + 1, row_length + 1);
    Q_SetEntry (qmatrix, row + 1, 0, row + 1, this->elements[row]);
    for (index = 0; index < row_length; ++index) {
      position = this->indices[row] + index;
      Q_SetEntry (qmatrix, row + 1, index + 1, this->indices[position] + 1, this->elements[position]);
    }
  }
}
#endif /* NO_LASPACK */

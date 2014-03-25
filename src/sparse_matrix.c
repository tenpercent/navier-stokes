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

  this->filled_row = -1;
  this->filled_element = size + 1;
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

  ++this->filled_element;
  while (this->filled_row < row) {
    ++(this->filled_row);
    this->indices[this->filled_row] = this->filled_element;
  }
  this->indices[this->filled_element] = col;
  this->elements[this->filled_element] = value;
}

void Sparse_matrix_Finish_filling (
    Sparse_matrix * this) {

  while (this->filled_row < this->size - 1) {
    ++(this->filled_row);
    this->indices[this->filled_row] = this->filled_element;
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

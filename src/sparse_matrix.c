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

#include "relax_system.h"

void relax_system (Sparse_matrix * lh_side, 
                   double * rh_side,
                   double relaxation_constant) {

  register unsigned i = 0;
  const unsigned dimension = lh_side->size;
  const unsigned total_nonzero_elements = lh_side->indices[dimension];

  for (i = 0; i < dimension; ++i) {
    rh_side[i] *= relaxation_constant;
  }

  for (i = 0; i < total_nonzero_elements; ++i) {
    lh_side->elements[i] *= relaxation_constant;
  }

  return;
}

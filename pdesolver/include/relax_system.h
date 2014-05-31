#include "sparse_matrix.h"


/* Multiply linear equations system's
 * left and right sides
 * by relaxation constant
 */
void relax_system (Sparse_matrix * lh_side, 
                   double * rh_side,
                   double relaxation_constant);

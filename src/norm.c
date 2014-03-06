#include "norm.h"
#include <math.h>

double residual_norm_C (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double)) {
  double norm = fabs (approximation[0] - approximated (time_upper_boundary, space_coordinate[0]));
  for (unsigned i = 1; i < dimension; ++i) {
    norm = max (fabs (approximation[i] - approximated (time_upper_boundary, space_coordinate[i])), norm);
  }
  return norm;
}

double residual_norm_L2 (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double)) {
  double norm = 0.;
  for (unsigned i = 0; i < dimension; ++i) {
    norm += square (approximation[i] - approximated (time_upper_boundary, space_coordinate[i]));
  }
  norm = sqrt (norm);
  return norm;
}

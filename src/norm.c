#include "norm.h"
#include <math.h>

inline double max (double x, double y) {
  return (x > y) ? x : y;
}

inline double square (double x) {
  return x * x;
}

double residual_norm_C (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double)) {
  double norm = fabs (approximation[0] - approximated (time_upper_boundary, space_coordinate[0]));
  unsigned space_step;
  for (space_step = 1; space_step < dimension; ++space_step) {
    norm = max (fabs (approximation[space_step] - approximated (time_upper_boundary, space_coordinate[space_step])), norm);
  }
  return norm;
}

double residual_norm_L2 (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double)) {
  double norm = 0.;
  unsigned space_step;
  for (space_step = 0; space_step < dimension; ++space_step) {
    norm += square (approximation[space_step] - approximated (time_upper_boundary, space_coordinate[space_step]));
  }
  norm = sqrt (norm);
  return norm;
}

#include "structures.h"

double residual_norm_C (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

double residual_norm_L2 (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

inline double max (double x, double y) {
  return x > y ? x : y;
}

inline double square (double x) {
  return x * x;
}
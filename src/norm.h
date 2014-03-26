// C[0, space_upper_boundary] norm of residual
// on time layer corresponding to
// time_upper_boundary
double residual_norm_C (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

// L2[0, space_upper_boundary] norm of residual
// on time layer corresponding to
// time_upper_boundary
double residual_norm_L2 (double * approximation, 
                        unsigned dimension, 
                        double * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

double max (double x, double y);
double square (double x);

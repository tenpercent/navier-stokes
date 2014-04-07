// C[0, space_upper_boundary] norm of residual
// on time layer corresponding to
// time_upper_boundary
double residual_norm_C (double const * approximation, 
                        unsigned dimension, 
                        double const * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

// L2[0, space_upper_boundary] norm of residual
// on time layer corresponding to
// time_upper_boundary
double residual_norm_L2 (double const * approximation, 
                        unsigned dimension, 
                        double const * space_coordinate, 
                        double time_upper_boundary,
                        double (*approximated) (double, double));

// maximum value of (x) and (y)
double max (double x, double y);

// (x) squared
double square (double x);

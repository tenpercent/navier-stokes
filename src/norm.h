#ifndef NORM
#define NORM

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

// find norm of ({function} ({node_values})) in C
// e. g. node_values = G(space_coordinates[space_step]),
// space_step in range (0, total_space_nodes);
// function(x) = exp(x)
double function_norm_C (double const * node_values,
                        unsigned dimension,
                        double (*function) (double));

// maximum value of (x) and (y)
double max (double x, double y);

// (x) squared
double square (double x);

#endif

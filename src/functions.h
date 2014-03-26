#ifndef FUNCTIONS
#define FUNCTIONS

#include <math.h>
#include "structures.h"

// u (x, t) from the pdf with problem description
double u_exact (double x, double t);

// rho (x, t) from the pdf with problem description
double rho_exact (double x, double t);

// g (x, t) from the pdf with problem description
double g_exact (double x, double t);

// (du / dt) (x, t)
double u_t (double x, double t);

// (du / dx) (x, t)
double u_x (double x, double t);

// (d(rho) / dt) (x, t)
double rho_t (double x, double t);

// (d(rho) / dx) (x, t)
double rho_x (double x, double t);

// ((d^2)u / (dx)^2) (x, t)
double u_xx (double x, double t);

// the result of application
// of the first equation of PDE system
// to (u_exact) and (g_exact)
double rhs_1st_equation (double x, double t, Gas_parameters * parameters);

// the result of application
// of the second equation of PDE system
// to (u_exact) and (g_exact)
double rhs_2nd_equation (double x, double t, Gas_parameters * parameters);

#endif

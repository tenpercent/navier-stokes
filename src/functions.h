#ifndef FUNCTIONS
#define FUNCTIONS

#include <math.h>
#include "structures.h"

double u_exact (double x, double t);

double ro_exact (double x, double t);

double log_ro (double x, double t);

double u_t (double x, double t);

double u_x (double x, double t);

double ro_t (double x, double t);

double ro_x (double x, double t);

double u_xx (double x, double t);

double rhs_1st_eqation (double x, double t, Gas_parameters * parameters);

double rhs_2nd_equation (double x, double t, Gas_parameters * parameters);

#endif

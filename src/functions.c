#ifndef FUNCTIONS
#define FUNCTIONS

#include <math.h>
#include "structures.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950288
#endif

double u_exact (double x, double t) {
  return cos(2 * M_PI * t) * sin(M_PI * x * x * .01);
}

double rho_exact (double x, double t) {
  return exp(t) * (cos (M_PI * x * .1) + 1.5);
}

double g_exact (double x, double t) {
  return t + log (cos (M_PI * x * .1) + 1.5);
}

double g_x (double x, double t) {
  return -.1 * M_PI * sin (.1 * M_PI * x) / (cos (.1 * M_PI * x) + 1.5);
}

double g_t (double x, double t) {
  return 1.;
}

double u_t (double x, double t) {
  return -2. * M_PI * sin(2 * M_PI * t) * sin(M_PI * x * x * .01);;
}

double u_x (double x, double t) {
  return .02 * M_PI * x * cos(2 * M_PI * t) * cos(M_PI * x * x * .01);
}

double rho_t (double x, double t) {
  return ro_exact (x, t);
}

double rho_x (double x, double t) {
  return -.1 * M_PI * exp(t) * sin (M_PI * x * .1);
}

double u_xx (double x, double t) {
  return .02 * M_PI * cos(2 * M_PI * t) * 
    ( cos(M_PI * x * x * .01) + 
      M_PI * -.02 * x * x * sin(M_PI * x * x * .01)
    );
}

double rhs_1st_equation (double x, double t, Gas_parameters * parameters) {
  return g_t (x, t) + 
          .5 * (
              2 * u_exact (x, t) * g_x (x, t) + 
              g_exact(x, t) * u_x (x, t) + 
              (2 - g_exact (x, t)) * u_x (x, t));
}

double rhs_2nd_equation (double x, double t, Gas_parameters * parameters) {
  return u_t (x, t) + 
          u_exact (x, t) * u_x (x, t) + 
          parameters->p_ro * g_x (x, t) -
          parameters->viscosity * exp (-g_exact (x, t)) * u_xx (x, t);
}

#endif

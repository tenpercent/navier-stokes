#include <math.h>
#include <structures.h>

double u_exact (double x, double t) {
  return cos(2 * M_PI * t) * sin(M_PI * x * x * .01);
}

double ro_exact (double x, double t) {
  return exp(t) * (cos (M_PI * x * .1) + 1.5);
}

double u_t (double x, double t) {
  return -2. * M_PI * sin(2 * M_PI * t) * sin(M_PI * x * x * .01);;
}

double u_x (double x, double t) {
  return .02 * M_PI * x * cos(2 * M_PI * t) * cos(M_PI * x * x * .01);
}

double ro_t (double x, double t) {
  return ro_exact (x, t);
}

double ro_x (double x, double t) {
  return -.1 * M_PI * exp(t) * sin (M_PI * x * .1);
}

double u_xx (double x, double t) {
  return .02 * M_PI * cos(2 * M_PI * t) * 
    ( cos(M_PI * x * x * .01) + 
      M_PI * -.02 * x * x * sin(M_PI * x * x * .01)
    );
}

double rhs_1st_eqation (double x, double t, Gas_parameters *) {
  return ro_t (x, t) + 
    ro_x (x, t) * u_exact (x, t) + 
    ro_exact (x, t) * u_x (x, t);
}

double rhs_2nd_equation (double x, double t, Gas_parameters * parameters) {
  return ro_exact (x, t) * u_t (x, t) + 
    ro_exact (x, t) * u_exact (x, t) * u_x (x, t) +
    parameters->p_ro - 
    parameters->viscosity * u_xx (x, t);
}

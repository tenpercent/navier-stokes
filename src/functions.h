#include <math.h>

// Dunno why. Really

double u_exact (double x, double t) {
  return cos (2 * M_PI * t) * sin (M_PI * x * x / 100);
}

double logp_exact (double x, double t) {
  return t + log (cos (M_PI * x / 10) + 1.5);
}

double v_rhs (double x, double t) {
  return 0.;
}

double g_rhs (double x, double t) {
  return 0.;
}

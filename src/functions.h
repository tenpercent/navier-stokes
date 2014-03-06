#include <math.h>

// Dunno why. Really

double u_bound (double x, double t) {
  return cos (2 * M_PI * t) * sin (M_PI * x * x / 100);
}

double logp_bound (double x, double t) {
  return t + log (cos (M_PI * x / 10) + 1.5);
}

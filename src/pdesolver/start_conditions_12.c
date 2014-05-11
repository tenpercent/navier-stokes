#include <math.h>
#include "start_conditions.h"

double u_start (double x) {
  return sin(M_PI * x * x * .01);
}

double rho_start (double x) {
  return cos (M_PI * x * .1) + 1.5;
}

double g_start (double x) {
  return log (cos (M_PI * x * .1) + 1.5);
}

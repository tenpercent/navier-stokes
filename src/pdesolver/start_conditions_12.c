#include <math.h>
#include "start_conditions.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950288
#endif

double u_start (double x) {
  return sin(M_PI * x * x * .01);
}

double rho_start (double x) {
  return cos (M_PI * x * .1) + 1.5;
}

double g_start (double x) {
  return log (cos (M_PI * x * .1) + 1.5);
}

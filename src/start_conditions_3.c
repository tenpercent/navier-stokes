#include <math.h>
#include <assert.h>
#include "start_conditions.h"

double u_start (double x) {

  double const epsilon = 1e-8;

  assert ((x > -epsilon) && (x < 10 + epsilon));

  return 0;
}

double rho_start (double x) {

  double const epsilon = 1e-8;

  assert ((x > -epsilon) && (x < 10 + epsilon));

  return (x > 4.5 - epsilon && x < 5.5 + epsilon) ? 2. : 1.;
}

double g_start (double x) {

  return log (rho_start (x));
}

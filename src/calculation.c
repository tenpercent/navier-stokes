#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "calculation.h"
#include "fill_system.h"
#include "start_conditions.h"
#include "functions.h"
#include "print.h"
#include "norm.h"

/* G and V values are packed as follows:
 *
 * G[0] V[0] G[1] V[1] ... G[M] V[M]
 *  1 .. 2 .. 3 .. 4 ..... 2M-1  2M
 */

#define G_INDEX(i) (2 * (i))
#define V_INDEX(i) (2 * (i) + 1)

/* Useful for arrays construction */
#define NEW(type, len) ((type*)malloc((len)*sizeof(type)))

#ifndef NO_LASPACK
#include <laspack/rtc.h>
void call_laspack_method (
    IterProcType    method,
    Sparse_matrix * lh_side,
    double        * unknown_vector,
    double        * rh_side,
    double          max_iterations,
    PrecondProcType preconditioner,
    double          relaxation_constant) {

  QMatrix lh_side_laspack;
  Vector rh_side_laspack;
  Vector unknown_vector_laspack;

  Q_Constr (&lh_side_laspack,
            "Matrix",
            lh_side->size, /* side */
            False, /* non-symmetric */
            Rowws, /* row-wise storage */
            Normal,
            True);
  V_Constr (&rh_side_laspack,
            "Right-hand side of equation",
            lh_side->size,
            Normal,
            False); /* we pass our own data */
  V_Constr (&unknown_vector_laspack,
            "Unknown vector",
            lh_side->size,
            Normal,
            False); /* the same here */

  Sparse_matrix_to_QMatrix (lh_side, &lh_side_laspack);

  /* Laspack Vector indices start with 1 */
  rh_side_laspack.Cmp = rh_side - 1;
  unknown_vector_laspack.Cmp = unknown_vector - 1;

  (*method) (&lh_side_laspack,
             &unknown_vector_laspack,
             &rh_side_laspack,
             max_iterations,
             preconditioner,
             relaxation_constant);

  Q_Destr (&lh_side_laspack);
  V_Destr (&rh_side_laspack);
  V_Destr (&unknown_vector_laspack);

}
#endif /* NO_LASPACK */

void find_approximate_solution (
    double * G,
    double * V,
    Node_status const * node_status,
    double const * space_coordinates,
    Gas_parameters const * gas_parameters,
    Grid const * grid,
    Iterative_Method_parameters const * iterative_method_parameters) {

  // left-hand side of the system
  Sparse_matrix lh_side;
  // right-hand side of the system
  double * rh_side = NEW(double, 2 * grid->X_nodes);
  // will be found when we solve the system
  double * unknown_vector = NEW(double, 2 * grid->X_nodes);
  // buffer used by iterative method
  double * buffer = NEW(double, 20 * grid->X_nodes);

  register unsigned time_step;

  unsigned const max_iterations = 50;
  Sparse_matrix_Construct (&lh_side, 2 * grid->X_nodes, 10 * grid->X_nodes - 10);

#ifndef NO_LASPACK
  // iterative algorithm accuracy
  SetRTCAccuracy (iterative_method_parameters->accuracy);
#endif /* NO_LASPACK */

  fill_unknown_vector (G, V, unknown_vector, grid->X_nodes);

  for (time_step = 1; time_step < grid->T_nodes; ++time_step) {
    // this may slow things up
    print_info_about_current_iteration (time_step, grid);

    fill_system (&lh_side, rh_side, grid, node_status, gas_parameters, space_coordinates, time_step, G, V);

    // launch iteration algorithm
    if (iterative_method_parameters->implementation == Implementation_Native) {
      iterative_method_parameters->method (&lh_side,
              unknown_vector,
              rh_side,
              max_iterations,
              iterative_method_parameters->preconditioner_type,
              iterative_method_parameters->relaxation_constant,
              iterative_method_parameters->accuracy,
              buffer);
    }
#ifndef NO_LASPACK
    else {
      call_laspack_method (
              iterative_method_parameters->method_laspack,
              &lh_side,
              unknown_vector,
              rh_side,
              max_iterations,
              iterative_method_parameters->preconditioner_type_laspack,
              iterative_method_parameters->relaxation_constant);
    }
#endif /* NO_LASPACK */

    fill_approximation (G, V, unknown_vector, grid->X_nodes);
  }

  Sparse_matrix_Destruct (&lh_side);
  free (unknown_vector);
  free (rh_side);
  free (buffer);

  return;
}

void fill_mesh_at_initial_time (
    double * G,
    double * V,
    double (*g) (double), // log (ro_0)
    double (*v) (double), // u_0
    double const * space_coordinates,
    unsigned space_nodes) {
  register unsigned space_step = 0;
  for (space_step = 0; space_step < space_nodes; ++space_step) {
    G[space_step] = g(space_coordinates[space_step]);
    V[space_step] = v(space_coordinates[space_step]);
  }
  return;
}

void fill_unknown_vector (
    double * G,
    double * V,
    double * unknown_vector,
    unsigned space_nodes) {

  register unsigned space_step = 0;
  for (space_step = 0; space_step < space_nodes; ++space_step) {
    unknown_vector[G_INDEX(space_step)] = G[space_step];
    unknown_vector[V_INDEX(space_step)] = V[space_step];
  }
  return;
}

void fill_approximation (double * G, double * V, double const * solutions, unsigned total_values) {
  register unsigned space_step = 0;

  for (space_step = 0; space_step < total_values; ++space_step) {
    G[space_step] = solutions[G_INDEX(space_step)];
    V[space_step] = solutions[V_INDEX(space_step)];
  }
  V[0] = V[total_values - 1] = 0.;

  return;
}

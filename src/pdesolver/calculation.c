#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "calculation.h"
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

/* Useful for Sparse_matrix filling */
#define MATRIX_APPEND(col, value) Sparse_matrix_Append_element (lh_side, row_number, col, value)

#ifndef NO_LASPACK
#include <rtc.h>
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
    Node_status * node_status, 
    double * space_coordinates, 
    Gas_parameters * gas_parameters,
    Grid * grid,
    Iterative_Method_parameters * iterative_method_parameters) {

  // left-hand side of the system
  Sparse_matrix lh_side; 
  // right-hand side of the system
  double * rh_side = NEW(double, 2 * grid->X_nodes); 
  // will be found when we solve the system
  double * unknown_vector = NEW(double, 2 * grid->X_nodes);
  // buffer used by iterative method
  double * buffer = NEW(double, 20 * grid->X_nodes);

  register unsigned time_step = 0,
           space_step = 0;

  unsigned max_iterations = 500;
  Sparse_matrix_Construct (&lh_side, 2 * grid->X_nodes, 10 * grid->X_nodes - 10);

  // initialize unknown vector
  // with the approximation of next time layer values
  // which is this time layer values
  // as the functions are continuous
  for (space_step = 0; space_step < grid->X_nodes; ++space_step) {
    unknown_vector[G_INDEX(space_step)] = g_exact (space_coordinates[space_step], 0.);
    unknown_vector[V_INDEX(space_step)] = u_exact (space_coordinates[space_step], 0.);
  }

#ifndef NO_LASPACK
  // iterative algorithm accuracy
  SetRTCAccuracy (iterative_method_parameters->accuracy);
#endif /* NO_LASPACK */

  fill_mesh_at_initial_time (G, V, g_exact, u_exact, space_coordinates, grid->X_nodes); 

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

void fill_system (
    Sparse_matrix * lh_side,
    double * rh_side,
    Grid const * grid,
    Node_status const * node_status,
    Gas_parameters * gas_parameters,
    double const * space_coordinates,
    unsigned const time_step,
    double const * G,
    double const * V) {

  register unsigned space_step = 0,
  // iterator through rows
           row_number = 0;

  // auxiliary constants
  double const h   = 1. / grid->X_step,
         h_2 = .5 * h,
         h_4 = .25 * h,
         h_6 = h / 6.,
         hh_4_3 = h * h * 4. / 3.,

         tau = 1. / grid->T_step;

  double const viscosity_norm = 
    gas_parameters->viscosity * function_norm_C (G, grid->X_nodes, exp_1);

  double const mu_tilda_hh_4_3 = viscosity_norm * hh_4_3; 

  Sparse_matrix_Clear (lh_side);

  for (space_step = 0; space_step < grid->X_nodes; ++space_step) {
    switch (node_status[space_step]) {
      case LEFT:
        MATRIX_APPEND (G_INDEX(0), tau - h_2 * V[0]);
        MATRIX_APPEND (V_INDEX(0), -h);
        MATRIX_APPEND (G_INDEX(1), h_2 * V[1]);
        MATRIX_APPEND (V_INDEX(1), h);

        rh_side[row_number] = tau * G[0] +
                   h_2 * G[0] * (V[1] - V[0]) +
                   h_4 *
                     (G[2] * V[2] - 2 * G[1] * V[1] + G[0] * V[0] +
                       (2 - G[0]) * (V[2] - 2 * V[1] + V[0])) +
                   rhs_1st_equation (space_coordinates[space_step], time_step * grid->T_step, gas_parameters);
        ++row_number;

        /* dummy equation */

        MATRIX_APPEND (V_INDEX(0), 1);
        rh_side[row_number] = 0;
        ++row_number;
        break;

      case MIDDLE:
        MATRIX_APPEND (G_INDEX(space_step - 1), -h_4 * V[space_step] - h_4 * V[space_step - 1]);
        MATRIX_APPEND (V_INDEX(space_step - 1), -h_2);
        MATRIX_APPEND (G_INDEX(space_step), tau);
        MATRIX_APPEND (G_INDEX(space_step + 1), h_4 * V[space_step] + h_4 * V[space_step + 1]);
        MATRIX_APPEND (V_INDEX(space_step + 1), h_2);

        rh_side[row_number] = tau * G[space_step] +
                   h_4 * G[space_step] * (V[space_step + 1] - V[space_step - 1]) +
                   rhs_1st_equation (space_coordinates[space_step], time_step * grid->T_step, gas_parameters);
        ++row_number;
        /* another equation */
        /* attention: these values should be changed when (viscosity != 0) */
        MATRIX_APPEND (G_INDEX(space_step - 1), -h_2 * gas_parameters->p_ro);
        // change me!
        MATRIX_APPEND (V_INDEX(space_step - 1), -h_6 * V[space_step] - h_6 * V[space_step - 1] - mu_tilda_hh_4_3);
        // change me!
        MATRIX_APPEND (V_INDEX(space_step), tau + 2 * mu_tilda_hh_4_3);
        MATRIX_APPEND (G_INDEX(space_step + 1), h_2 * gas_parameters->p_ro);
        // change me!
        MATRIX_APPEND (V_INDEX(space_step + 1), h_6 * V[space_step] + h_6 * V[space_step + 1] - mu_tilda_hh_4_3);
        // change me!
        rh_side[row_number] = tau * V[space_step] -
                  hh_4_3 * 
                    (viscosity_norm - gas_parameters->viscosity * exp_1 (G[space_step])) *
                    (V[space_step - 1] - 2 * V[space_step] + V[space_step + 1]) + 
                  rhs_2nd_equation (space_coordinates[space_step], time_step * grid->T_step, gas_parameters);
        ++row_number;
        break;

      case RIGHT:
        MATRIX_APPEND (G_INDEX(grid->X_nodes - 2), -h_2 * V[grid->X_nodes - 2]);
        MATRIX_APPEND (V_INDEX(grid->X_nodes - 2), -h);
        MATRIX_APPEND (G_INDEX(grid->X_nodes - 1), tau + h_2 * V[grid->X_nodes - 1]);
        MATRIX_APPEND (V_INDEX(grid->X_nodes - 1), h);

        rh_side[row_number] = tau * G[grid->X_nodes - 1] +
                   h_2 * G[grid->X_nodes - 1] * (V[grid->X_nodes - 1] - V[grid->X_nodes - 2]) -
                   h_4 * 
                     (G[grid->X_nodes - 1] * V[grid->X_nodes - 1] - 
                      2 * G[grid->X_nodes - 2] * V[grid->X_nodes - 2] + 
                      G[grid->X_nodes - 3] * V[grid->X_nodes - 3] +
                        (2 - G[grid->X_nodes - 1]) * 
                        (V[grid->X_nodes - 1] - 2 * V[grid->X_nodes - 2] + V[grid->X_nodes - 3])) +
                   rhs_1st_equation(space_coordinates[space_step], time_step * grid->T_step, gas_parameters);
        ++row_number;

        /* dummy equation */ 

        MATRIX_APPEND (V_INDEX(grid->X_nodes - 1), 1);
        rh_side[row_number] = 0;
        ++row_number;
        break;

      default:
        break;
    }
  }

  Sparse_matrix_Finish_filling (lh_side);
  return;
}

void fill_mesh_at_initial_time (
    double * G,
    double * V,
    double (*g) (double, double), // log (ro_0)
    double (*v) (double, double), // u_0
    double * space_coordinates,
    unsigned space_nodes) {
  register unsigned space_step = 0;
  for (space_step = 0; space_step < space_nodes; ++space_step) {
    G[space_step] = g(space_coordinates[space_step], 0);
    V[space_step] = v(space_coordinates[space_step], 0);
  }
  return;
}

void fill_approximation (double * G, double * V, double * solutions, unsigned total_values) {
  register unsigned space_step = 0;

  for (space_step = 0; space_step < total_values; ++space_step) {
    G[space_step] = solutions[G_INDEX(space_step)];
    V[space_step] = solutions[V_INDEX(space_step)];
  }
  V[0] = V[total_values - 1] = 0.;

  return;
}
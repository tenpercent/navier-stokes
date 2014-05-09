#include <math.h>
#include "fill_system.h"
#include "functions.h"
#include "norm.h"

/* G and V values are packed as follows:
 *
 * G[0] V[0] G[1] V[1] ... G[M] V[M]
 *  1 .. 2 .. 3 .. 4 ..... 2M-1  2M
 */
#define G_INDEX(i) (2 * (i))
#define V_INDEX(i) (2 * (i) + 1)

/* Useful for Sparse_matrix filling */
#define MATRIX_APPEND(col, value) Sparse_matrix_Append_element (lh_side, row_number, col, value)

void fill_system (
    Sparse_matrix * lh_side,
    double * rh_side,
    Grid const * grid,
    Node_status const * node_status,
    Gas_parameters const * gas_parameters,
    double const * space_coordinates,
    unsigned const time_step,
    double const * G,
    double const * V) {

  unsigned space_step = 0,
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

  // printf ("scaled norm of mutilda is %lf", mu_tilda_hh_4_3);
  // printf ("hh43 is %lf", hh_4_3);

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

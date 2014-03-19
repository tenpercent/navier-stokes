#include "calculation.h"
#include "functions.h"

/* G and V values are packed as follows:
 *
 * G[0] V[0] G[1] V[1] ... G[M] V[M]
 * 
 * I cannot into macros, so working fix would be welcome
 */

#define G_INDEX(i) (2 * i + 1)
#define V_INDEX(i) (2 * i + 2)

void next_TimeLayer_Calculate (
    double * G, 
    double * V, 
    Node_status * node_status, 
    double * space_coordinates, 
    Gas_parameters * parameters,
    Grid * grid) {
  QMatrix lh_side; // left-hand side of the system
  Vector rh_side; // right-hand side of the system
  Vector unknown_vector; // which will be found when we solve the system
  unsigned time_step;

  Q_Constr (&lh_side, 
            "Matrix", 
            2 * grid->X_nodes, // side
            False, // non-symmetric
            Rowws, // row-wise storage
            Normal, // internal laspack stuff
            True); // internal laspack stuff
  V_Constr (&rh_side, 
            "Right-hand side of equation", 
            2 * grid->X_nodes, // size
            Normal, // internal laspack stuff
            True); // internal laspack stuff
  V_Constr (&unknown_vector, 
            "Unknown vector", 
            2 * grid->X_nodes, 
            Normal, 
            True);
  SetRTCAccuracy (1e-8);

  // T == 0
  fill_mesh_at_initial_time (G, V, log_ro, u_exact, space_coordinates, grid->X_nodes); 

  for (time_step = 1; time_step < grid->T_nodes; ++time_step) {
    fill_system (&lh_side, &rh_side, grid, node_status, parameters, space_coordinates, time_step, G, V);

    // launch iteration algorithm
    CGNIter (&lh_side, 
            &unknown_vector, 
            &rh_side, 
            2000, // max iterations
            JacobiPrecond, // preconditioner type
            1.2); // preconditioner relaxation constant; probably, should be changed

    fill_approximation (G, V, &unknown_vector);
  }

  Q_Destr (&lh_side);
  V_Destr (&unknown_vector);
  V_Destr (&rh_side);

  return;
}

void set_qmatrix_entries (
    QMatrix * matrix, 
    unsigned row, 
    unsigned * nonzero_columns, 
    double * values, 
    unsigned row_length) {
  unsigned entry_number;
  for (entry_number = 0; entry_number < row_length; ++entry_number) {
    Q_SetEntry (matrix, row, entry_number, nonzero_columns[entry_number], values[entry_number]);
  }
  return;
}

void fill_system (
    QMatrix * lh_side,
    Vector * rh_side,
    Grid * grid,
    Node_status * node_status,
    Gas_parameters * parameters,
    double * space_coordinates,
    unsigned time_step,
    double * G,
    double * V) {
  unsigned space_step = 0;
  unsigned first_type_equation_coef_length = 5;
  unsigned second_type_equation_coef_length = 4;
  unsigned third_type_equation_coef_length = 4;
  unsigned fourth_type_equation_coef_length = 5;
  unsigned equation_number = 1;

  unsigned nonzero_columns[5];
  double lh_values[5];
  double rh_value;

  for (space_step = 0; space_step < grid->X_nodes; ++space_step) {
    switch (node_status[space_step]) {
      case LEFT:
        nonzero_columns[0] = G_INDEX(0); // G_0
        nonzero_columns[1] = V_INDEX(0); // V_0
        nonzero_columns[2] = G_INDEX(1); // G_1
        nonzero_columns[3] = V_INDEX(1); // V_1

        lh_values[0] = grid->T_step - .5 * grid->X_step * V[0];
        lh_values[1] = -grid->X_step;
        lh_values[2] = .5 * grid->X_step * V[1];
        lh_values[3] = grid->X_step;

        rh_value = grid->T_step * G[0] +
                   grid->X_step * G[0] * (V[1] - V[0]) +
                   .25 * grid->X_step *
                     (G[2] * V[2] - 2 * G[1] * V[1] + G[0] * V[0] +
                      (2 - G[0]) * (V[2] - 2 * V[1] + V[0])) +
                   rhs_1st_equation (space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, second_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, second_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;
        break;

      case MIDDLE:
        nonzero_columns[0] = G_INDEX(space_step - 1); // G_{n-1}
        nonzero_columns[1] = V_INDEX(space_step - 1); // V_{n-1}
        nonzero_columns[2] = G_INDEX(space_step);     // G_{n}
        nonzero_columns[3] = G_INDEX(space_step + 1); // G_{n+1}
        nonzero_columns[4] = V_INDEX(space_step + 1); // V_{n+1}

        lh_values[0] = -.25 * grid->X_step * V[space_step] - .25 * grid->X_step * V[space_step - 1];
        lh_values[1] = -grid->X_step;
        lh_values[2] = grid->T_step;
        lh_values[3] = .25 * grid->X_step * V[space_step] + .25 * grid->X_step * V[space_step + 1];
        lh_values[4] = grid->X_step;

        rh_value = grid->T_step * G[space_step] +
                   .5 * grid->X_step * G[space_step] * (V[space_step + 1] - V[space_step - 1]) +
                   rhs_1st_equation (space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, first_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, first_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;

        /* another equation */

        nonzero_columns[0] = G_INDEX(space_step - 1); // G_{n-1}
        nonzero_columns[1] = V_INDEX(space_step - 1); // V_{n-1}
        nonzero_columns[2] = V_INDEX(space_step);     // V_{n}
        nonzero_columns[3] = G_INDEX(space_step + 1); // G_{n+1}
        nonzero_columns[4] = V_INDEX(space_step + 1); // V_{n+1}

        // attention: these values should be changed when (viscosity != 0)
        lh_values[0] = -.5 * grid->X_step * parameters->p_ro;
        lh_values[1] = -1. / 6 * grid->X_step * V[space_step] - 
                        1. / 6 * grid->X_step * V[space_step - 1];
        lh_values[2] = grid->T_step;
        lh_values[3] = .5 * grid->X_step * parameters->p_ro;
        lh_values[4] = 1. / 6 * grid->X_step * V[space_step] + 
                       1. / 6 * grid->X_step * V[space_step + 1];
        // and these
        rh_value = grid->T_step * V[space_step] +
                   rhs_2nd_equation (space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, fourth_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, fourth_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;
        break;

      case RIGHT:
        nonzero_columns[0] = G_INDEX(grid->X_nodes - 2); // G_{M-1}
        nonzero_columns[1] = V_INDEX(grid->X_nodes - 2); // V_{M-1}
        nonzero_columns[2] = G_INDEX(grid->X_nodes - 1); // G_{M}
        nonzero_columns[3] = V_INDEX(grid->X_nodes - 1); // V_{M}

        lh_values[0] = -.5 * grid->X_step * V[grid->X_nodes - 2];
        lh_values[1] = -grid->X_step;
        lh_values[2] = grid->T_step + .5 * grid->X_step * V[grid->X_nodes - 1];
        lh_values[3] = grid->X_step;

        rh_value = G[grid->X_nodes - 1] * grid->T_step +
                   .5 * grid->T_step * G[grid->X_nodes - 1] * (V[grid->X_nodes - 1] - V[grid->X_nodes - 2]) +
                   -.25 * grid->X_step * 
                     (G[grid->X_nodes - 1] * V[grid->X_nodes - 1] - 2 * G[grid->X_nodes - 2] * V[grid->X_nodes - 2] + G[grid->X_nodes - 3] * V[grid->X_nodes - 3] +
                        (2 - G[grid->X_nodes - 1]) * (V[grid->X_nodes - 1] - 2 * V[grid->X_nodes - 2] + V[grid->X_nodes - 3])) +
                   rhs_1st_equation(space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, third_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, third_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;
        break;

      default:
        break;
    }
  }
}

void fill_mesh_at_initial_time (
    double * G,
    double * V,
    double (*g) (double, double), // log (ro_0)
    double (*v) (double, double), // u_0
    double * space_coordinates,
    unsigned space_nodes) {
  unsigned space_step = 0;
  for (space_step = 0; space_step < space_nodes; ++space_step) {
    G[space_step] = g(space_coordinates[space_step], 0);
    V[space_step] = v(space_coordinates[space_step], 0);
  }
  return;
}

void fill_approximation (double * G, double * V, Vector * unknown_vector) {
  unsigned total_values = V_GetDim (unknown_vector) / 2;
  unsigned space_step = 0;
  for (space_step = 0; space_step < total_values; ++space_step) {
      G[space_step] = V_GetCmp (unknown_vector, G_INDEX(space_step));
    }
    V[0] = 0.;
    for (space_step = 1; space_step < total_values - 1; ++space_step) {
      V[space_step] = V_GetCmp (unknown_vector, V_INDEX(space_step));
    }
    V[total_values - 1] = 0.;
    return;
}

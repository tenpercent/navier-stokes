#include "calculation.h"
#include "functions.h"

/* G and V values are packed as follows:
 *
 * G[0] V[0] G[1] V[1] ... G[M] V[M]
 *  1 .. 2 .. 3 .. 4 ..... 2M-1  2M
 */

#define G_INDEX(i) (2 * (i) + 1)
#define V_INDEX(i) (2 * (i) + 2)

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
  fill_mesh_at_initial_time (G, V, g_exact, u_exact, space_coordinates, grid->X_nodes); 

  for (time_step = 1; time_step < grid->T_nodes; ++time_step) {
    printf ("\rTime step is %d of %d.", time_step, grid->T_nodes - 1);
    fflush (stdout);
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

  double _h = 1. / grid->X_step;
  double _h_2 = .5 * _h;
  double _h_4 = .25 * _h;
  double _h_6 = _h / 6.;

  double _tau = 1. / grid->T_step;

  for (space_step = 0; space_step < grid->X_nodes; ++space_step) {
    switch (node_status[space_step]) {
      case LEFT:
        nonzero_columns[0] = G_INDEX(0); // G_0
        nonzero_columns[1] = V_INDEX(0); // V_0
        nonzero_columns[2] = G_INDEX(1); // G_1
        nonzero_columns[3] = V_INDEX(1); // V_1

        lh_values[0] = _tau - _h_2 * V[0];
        lh_values[1] = -_h;
        lh_values[2] = _h_2 * V[1];
        lh_values[3] = _h;

        rh_value = _tau * G[0] +
                   _h_2 * G[0] * (V[1] - V[0]) +
                   _h_4 *
                     (G[2] * V[2] - 2 * G[1] * V[1] + G[0] * V[0] +
                      (2 - G[0]) * (V[2] - 2 * V[1] + V[0])) +
                   rhs_1st_equation (space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, second_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, second_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;

        /* dummy equation */

        nonzero_columns[0] = V_INDEX(0);
        lh_values[0] = 1;
        rh_value = 0;
        Q_SetLen (lh_side, equation_number, 1);
        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, 1);
        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;
        break;

      case MIDDLE:
        nonzero_columns[0] = G_INDEX(space_step - 1); // G_{n-1}
        nonzero_columns[1] = V_INDEX(space_step - 1); // V_{n-1}
        nonzero_columns[2] = G_INDEX(space_step);     // G_{n}
        nonzero_columns[3] = G_INDEX(space_step + 1); // G_{n+1}
        nonzero_columns[4] = V_INDEX(space_step + 1); // V_{n+1}

        lh_values[0] = -_h_4 * V[space_step] - _h_4 * V[space_step - 1];
        lh_values[1] = -_h;
        lh_values[2] = _tau;
        lh_values[3] = _h_4 * V[space_step] + _h_4 * V[space_step + 1];
        lh_values[4] = _h;

        rh_value = _tau * G[space_step] +
                   _h_2 * G[space_step] * (V[space_step + 1] - V[space_step - 1]) +
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
        lh_values[0] = -_h_2 * parameters->p_ro;
        lh_values[1] = -_h_6 * V[space_step] - 
                        _h_6 * V[space_step - 1];
        lh_values[2] = _tau;
        lh_values[3] = _h_2 * parameters->p_ro;
        lh_values[4] = _h_6 * V[space_step] + 
                       _h_6 * V[space_step + 1];
        // and these
        rh_value = _tau * V[space_step] +
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

        lh_values[0] = -_h_2 * V[grid->X_nodes - 2];
        lh_values[1] = -_h;
        lh_values[2] = _tau + _h_2 * V[grid->X_nodes - 1];
        lh_values[3] = _h;

        rh_value = _tau * G[grid->X_nodes - 1] +
                   _h_2 * G[grid->X_nodes - 1] * (V[grid->X_nodes - 1] - V[grid->X_nodes - 2]) -
                   _h_4 * 
                     (G[grid->X_nodes - 1] * V[grid->X_nodes - 1] - 2 * G[grid->X_nodes - 2] * V[grid->X_nodes - 2] + G[grid->X_nodes - 3] * V[grid->X_nodes - 3] +
                        (2 - G[grid->X_nodes - 1]) * (V[grid->X_nodes - 1] - 2 * V[grid->X_nodes - 2] + V[grid->X_nodes - 3])) +
                   rhs_1st_equation(space_coordinates[space_step], time_step * grid->T_step, parameters);

        Q_SetLen (lh_side, equation_number, third_type_equation_coef_length);

        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, third_type_equation_coef_length);

        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;

        /* dummy equation */ 
         
        nonzero_columns[0] = V_INDEX(grid->X_nodes - 1);
        lh_values[0] = 1;
        rh_value = 0;
        Q_SetLen (lh_side, equation_number, 1);
        set_qmatrix_entries (lh_side, equation_number, nonzero_columns, lh_values, 1);
        V_SetCmp (rh_side, equation_number, rh_value);
        ++equation_number;
        break;

      default:
        break;
    }
  }
  return;
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
  size_t total_values = V_GetDim (unknown_vector) >> 1;
  unsigned space_step = 0;

  for (space_step = 0; space_step < total_values; ++space_step) {
      G[space_step] = V_GetCmp (unknown_vector, G_INDEX(space_step));
  }

  V[0] = 0.;
  for (space_step = 1; space_step + 1 < total_values; ++space_step) { // total_values == 0 wtf
    V[space_step] = V_GetCmp (unknown_vector, V_INDEX(space_step));
  }
  V[total_values - 1] = 0.;
  return;
}

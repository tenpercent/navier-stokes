#include "calculation.h"

/* G and V values are packed as follows:
 *
 * G[0] G[1] ... G[M] V[1] V[2] ... V[M-1]
 *  0    1        M   M+1  M+2       2M-1
 *
 * These macros make it possible to easily change packing order in
 * the future.
 */
 
#define G_IND(i) (i)
#define V_IND(i) (grid->X_nodes - 1 + (i))

void next_TimeLayer_Calculate (
    double * G, 
    double * V, 
    Node_status * node_status, 
    double * space_coordinate, 
    Gas_parameters * parameters,
    Grid * grid) {
  QMatrix lh_side; // left-hand side of the system
  Vector rh_side; // right-hand side of the system
  Vector unknown_vector; // which will be found when we solve the system
  unsigned time_step, i;

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
  // fill_mesh_at_initial_time (...); 

  for (time_step = 1; time_step < grid->T_nodes; ++time_step) {
    fill_system (&lh_side, &rh_side, grid, node_status, parameters, G, V);

    // launch iteration algorithm
    CGNIter (&lh_side, 
            &unknown_vector, 
            &rh_side, 
            2000, // max iterations
            JacobiPrecond, // preconditioner type
            1.2); // preconditioner relaxation constant; probably, should be changed

    for (i = 0; i < grid->X_nodes; ++i) {
      G[i] = V_GetCmp (&unknown_vector, G_IND(i));
    }
    for (i = 1; i < grid->X_nodes - 1; ++i) {
      V[i] = V_GetCmp (&unknown_vector, V_IND(i));
    }
  }

  Q_Destr (&lh_side);
  V_Destr (&unknown_vector);
  V_Destr (&rh_side);

  return;
}

void fill_system (
    QMatrix * lh_side,
    Vector * rh_side,
    Grid * grid,
    Node_status * node_status,
    Gas_parameters * parameters,
    double * G,
    double * V) {

  unsigned m, M = grid->X_nodes - 1;
  double h = grid->X_step, tau = grid->T_step;
  int i = -1;

  Q_SetLen  (lh_side, G_IND(0), 4);
  Q_SetEntry(lh_side, G_IND(0), 0, G_IND(0), 1. / tau - .5 * V[0] / h);
  Q_SetEntry(lh_side, G_IND(0), 1, G_IND(1), .5 * V[1] / h);
  Q_SetEntry(lh_side, G_IND(0), 2, V_IND(0), -1. / h);
  Q_SetEntry(lh_side, G_IND(0), 3, V_IND(1), 1. / h);
  V_SetCmp  (rh_side, G_IND(0), G[0] / tau + G[0] * .5 * (V[1] - V[0]) / h +
                                .25 * (G[2] * V[2] - 2 * G[1] * V[1] + G[0] * V[0]) / h +
                                .25 * (V[2] - 2 * V[1] + V[0]) * (2 - G[0]) / h);

  for (m = 1; m < M; ++m) {
    Q_SetLen  (lh_side, G_IND(m), 5);

    Q_SetEntry(lh_side, G_IND(m), 0, G_IND(m - 1), -(V[m] + V[m-1]) * .25 / h);
    Q_SetEntry(lh_side, G_IND(m), 1, G_IND(m), G[m] / tau);
    Q_SetEntry(lh_side, G_IND(m), 2, G_IND(m + 1), (V[m] + V[m+1]) * .25 / h);
    Q_SetEntry(lh_side, G_IND(m), 3, V_IND(m - 1), -.5 / h);
    Q_SetEntry(lh_side, G_IND(m), 4, V_IND(m + 1), .5 / h);

    V_SetCmp  (rh_side, G_IND(m), G[m] / tau + G[m] * (V[m+1] - V[m-1]) * .25 / h);
  }

  Q_SetLen  (lh_side, G_IND(M), 4);

  Q_SetEntry(lh_side, G_IND(M), 0, G_IND(M - 1), -.5 * V[M-1] / h);
  Q_SetEntry(lh_side, G_IND(M), 1, G_IND(M), 1. / tau + .5 * V[M] / h);
  Q_SetEntry(lh_side, G_IND(M), 2, V_IND(M - 1), -1. / h);
  Q_SetEntry(lh_side, G_IND(M), 3, V_IND(M), 1. / h);

  V_SetCmp  (rh_side, G_IND(M), G[M] / tau + G[M] * .5 * (V[M] - V[M-1]) / h -
                                .25 * (G[M] * V[M] - 2 * G[M-1] * V[M-1] + G[M-2] * V[M-2]) / h -
                                .25 * (V[M] - 2 * V[M-1] + V[M-2]) * (2 - G[M]) / h);

  i = -1;
  for (m = 1; m < M; ++m) {
    Q_SetLen  (lh_side, V_IND(m), 3 + (m > 1) + (m < M - 1));

    Q_SetEntry(lh_side, V_IND(m), ++i, G_IND(m - 1), -.5 * parameters->p_ro / h);
    Q_SetEntry(lh_side, V_IND(m), ++i, G_IND(m + 1), .5 * parameters->p_ro / h);
    if (m > 1) {
      Q_SetEntry(lh_side, V_IND(m), ++i, V_IND(m - 1), -(V[m] + V[m-1]) / (h * 6));
    }
    Q_SetEntry(lh_side, V_IND(m), ++i, V_IND(m), 1. / tau);
    if (m < M - 1) {
      Q_SetEntry(lh_side, V_IND(m), ++i, V_IND(m + 1), (V[m] + V[m+1]) / (h * 6));
    }

    V_SetCmp  (rh_side, V_IND(m), V[m] / tau);
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

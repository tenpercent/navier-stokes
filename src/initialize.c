#include "initialize.h"
#include <string.h>

void gas_parameters_Initialize (Gas_parameters * parameters) {
  parameters->time_upper_boundary   = 1.;
  parameters->space_upper_boundary  = 10.;
  parameters->p_ro                  = 10.;
  parameters->viscosity             = 0.;

  return;
}

void grid_Initialize (
    Grid * grid, 
    Gas_parameters * parameters,
    unsigned time_steps_magnifier, 
    unsigned space_steps_magnifier) {
  grid->X_nodes = 400 << space_steps_magnifier;
  grid->T_nodes = 400 << time_steps_magnifier;

  grid->X_step = parameters->space_upper_boundary / (grid->X_nodes - 1);
  grid->T_step = parameters->time_upper_boundary / (grid->T_nodes - 1);

  return;
}

void mesh_Initialize (Node_status * node_status, double * space_coordinate, Grid * grid) {
  unsigned space_step;
  node_status[0] = LEFT;
  node_status[grid->X_nodes - 1] = RIGHT;
  for (space_step = 1; space_step < grid->X_nodes - 1; ++space_step) {
    node_status[space_step] = MIDDLE;
  }

  for (space_step = 0; space_step < grid->X_nodes; ++space_step) {
    space_coordinate[space_step] = space_step * grid->X_step;
  }

  return;
}

void initialize_iterative_algorithm_parameters (Preconditioner_type * type, Iterative_method * method, int argc, char ** argv) {
  if (argc < 3) {
    *type = Jacobi;
    *method = CGN;
    return;
  } else {
    if (strncmp (argv[1], "Jacobi", 6) == 0) {
      *type = Jacobi;
    } else if (strncmp (argv[1], "SSOR", 4) == 0) {
      *type = SSOR;
    } else {
      *type = Jacobi;
    }

    if (strncmp (argv[2], "CGN", 3) == 0) {
      *method = CGN;
    } else if (strncmp (argv[2], "BiCGStab", 8) == 0) {
      *method = BiCGStab;
    } else if (strncmp (argv[2], "CGS", 3) == 0) {
      *method = CGS;
    } else if (strncmp (argv[2], "QMR", 3) == 0) {
      *method = QMR;
    } else if (strncmp (argv[2], "GMRES", 3) == 0) {
      *method = GMRES;
    } else {
      *method = CGN;
    }

    return;
  }
}

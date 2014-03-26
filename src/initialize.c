#include "initialize.h"
#include <strings.h>
#include <itersolv.h>

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

void initialize_iterative_algorithm_parameters (Iterative_Method_parameters * parameters, int argc, char ** argv) {
   parameters->preconditioner_type = &Precond_Jacobi;
   parameters->method = &Iterative_method_BiCGSTAB;
/*
  if (argc < 3) {
    parameters->preconditioner_type = JacobiPrecond;
    parameters->method = CGNIter;
    return;
  } {
    if (strncasecmp (argv[1], "Jacobi", 6) == 0) {
      parameters->preconditioner_type = JacobiPrecond;
    } else if (strncasecmp (argv[1], "SSOR", 4) == 0) {
      parameters->preconditioner_type = SSORPrecond;
    } else {
      // default value
      parameters->preconditioner_type = JacobiPrecond;
    }

    if (strncasecmp (argv[2], "CGN", 3) == 0) {
      parameters->method = CGNIter;
    } else if (strncasecmp (argv[2], "BiCGStab", 8) == 0) {
      parameters->method = BiCGSTABIter;
    } else if (strncasecmp (argv[2], "CGS", 3) == 0) {
      parameters->method = CGSIter;
    } else if (strncasecmp (argv[2], "QMR", 3) == 0) {
      parameters->method = QMRIter;
    } else if (strncasecmp (argv[2], "GMRES", 5) == 0) {
      parameters->method = GMRESIter;
    } else {
      // default value
      parameters->method = CGNIter;
    }

    return;
  }
*/
}

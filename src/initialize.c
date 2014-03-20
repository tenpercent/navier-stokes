#include "initialize.h"

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
  grid->X_nodes = 100 << space_steps_magnifier;
  grid->T_nodes = 100 << time_steps_magnifier;

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

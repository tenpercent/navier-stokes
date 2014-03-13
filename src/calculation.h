#include "structures.h"

#include <vector.h>
#include <errhandl.h>
#include <qmatrix.h>
#include <itersolv.h>
#include <rtc.h>
#include <operats.h>

void next_TimeLayer_Calculate (
  double * G, 
  double * V, 
  Node_status * node_status, 
  double * space_coordinate, 
  Gas_parameters * parameters,
  Grid * grid);

void fill_system (
  QMatrix * lh_side,
  Vector * rh_side,
  Grid * grid,
  Node_status * node_status,
  Gas_parameters * parameters,
  double * G,
  double * V);

void fill_mesh_at_initial_time (
  double * G,
  double * V,
  double (*g) (double, double),
  double * space_coordinates,
  unsigned space_nodes);

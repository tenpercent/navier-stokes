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
  Grid * grid);

void fill_system (
    QMatrix * lh_side,
    Vector * rh_side,
    Grid * grid,
    Node_status * node_status,
    double * G,
    double * V);

#include "calculation.h"
#include <vector.h>
#include <errhandl.h>
#include <qmatrix.h>
#include <itersolv.h>
#include <rtc.h>
#include <operats.h>
#include <version.h>
#include <copyrght.h>

void next_time_layer_Calculate (
    double * G, 
    double * V, 
    Node_status * node_status, 
    double * space_coordinate, 
    Grid * grid) {
  QMatrix lh_side; // left-hand side of the system
  Vector rh_side; // right-hand side of the system
  Vector unknown_vector; // which will be found when we solve the system

  Q_Constr (&lh_side, 
            "Matrix", 
            2 * grid->X_nodes, // side
            False, // non-symmetric
            Rowws, // row-wise storage
            Normal, // internal laspack stuff
            True); // internal laspack stuff
  V_Constr (&rh_side, 
            "Right-hand side", 
            2 * grid->X_nodes, // size
            Normal, // internal laspack stuff
            True); // internal laspack stuff
  V_Constr (&unknown_vector, 
            "Unknown vector", 
            2 * grid->X_nodes, 
            Normal, 
            True);

  SetRTCAccuracy (1e-8);

  /*
  заполняем матрицу и правую часть
  */

  // launch iteration algorithm
  CGNIter (&lh_side, 
          &unknown_vector, 
          &rh_side, 
          2000, // max iterations
          JacobiPrecond, // preconditioner type
          1.2); // preconditioner relaxation constant; probably, should be changed
  return;
}

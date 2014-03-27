#include <stdlib.h>
#include "construction.h"

void results_Construct (double ** results, unsigned max_global_iteration) {
  unsigned i;
  for (i = 0; i < RESULTS_SIZE; ++i) {
    results[i] = (double *) malloc (max_global_iteration * sizeof (double));
  }

  return;
}

void results_Destruct (double ** results) {
  unsigned i;
  for (i = 0; i < RESULTS_SIZE; ++i) {
    free (results[i]);
  }
  return;
}

void scheme_elements_Construct (double ** G, double ** V, unsigned X_nodes) {
  *V = (double *) malloc (X_nodes * sizeof (double));
  *G = (double *) malloc (X_nodes * sizeof (double));
  return;
}

void scheme_elements_Destruct (double * G, double * V) {
  free (V);
  free (G);
  return;
}

void mesh_elements_Construct (Node_status ** node_statuses, double ** space_coordinates, unsigned X_nodes) {
  *space_coordinates = (double *) malloc (X_nodes * sizeof (double));
  *node_statuses = (Node_status *) malloc (X_nodes * sizeof (Node_status));
  return;
}

void mesh_elements_Destruct (Node_status * node_statuses, double * space_coordinates) {
  free (space_coordinates);
  free (node_statuses);
  return;
}

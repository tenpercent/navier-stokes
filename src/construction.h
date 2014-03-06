#include "structures.h"
#define RESULTS_SIZE 5

void results_Construct (double ** results, unsigned max_global_iteration);

void results_Destruct (double ** results);

void scheme_elements_Construct (double * G, double * V, unsigned X_nodes);

void scheme_elements_Destruct (double * G, double * V);

void mesh_elements_Construct (Node_status * node_statuses, double * space_coordinates, unsigned X_nodes);

void mesh_elements_Destruct (Node_status * node_statuses, double * space_coordinates);
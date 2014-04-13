#ifndef CONSTRUCTION
#define CONSTRUCTION

#include "structures.h"
#define RESULTS_SIZE 8

// allocate memory for results
void results_Construct (double ** results, unsigned max_global_iteration);

// free memory of results
void results_Destruct (double ** results);

// allocate memory for arrays V and G
// which store the values of sought functions
// in mesh nodes
void scheme_elements_Construct (double ** G, double ** V, unsigned X_nodes);

// free memory of arrays V and G
void scheme_elements_Destruct (double * G, double * V);

// allocate memory for the information of node positions in the mesh (node_statuses)
// and the coordinates of the mesh nodes (space_coordinates)
void mesh_elements_Construct (Node_status ** node_statuses, double ** space_coordinates, unsigned X_nodes);

// free memory of the information of node positions in the mesh
// and the coordinates of the mesh nodes
void mesh_elements_Destruct (Node_status * node_statuses, double * space_coordinates);

#endif

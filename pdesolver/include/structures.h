#ifndef STRUCTURES
#define STRUCTURES

/* sctructure to store the problem parameters */
typedef struct {
  double time_upper_boundary;     /* time boundaries: [0, T]                                          */
  double space_upper_boundary;    /* space boundaries: [0, X]                                         */
  double p_ro;                    /* p(ro) is the function from the problem statement                 */
                                  /* we consider it is linear                                         */
                                  /* i.e. p(ro) = p_ro * ro                                           */
  double viscosity;               /* viscosity; also look into the problem statement. MU              */
  double artificial_viscosity;    /* artificial viscosity; also look into the problem statement. eta  */
} Gas_parameters;


/* structure to store grid parameters: 
 * nodes count and 
 * distance between neighbor nodes
 * in time and space directions
 */
typedef struct {
  unsigned X_nodes; /* total space nodes; should be >= 3. */
  unsigned T_nodes; /* total time nodes                   */
  double X_step;    /* distance between space nodes       */
  double T_step;    /* distance between time nodes        */
} Grid;


/* equations at the edges of space segment and in its middle are different. */
/* because of chosen differential scheme                                    */
/* so we need to know where we are when we construct the matrix of SLE      */
/* corresponding to chosen differential scheme                              */
typedef enum {
  LEFT = 0,
  MIDDLE = 1,
  RIGHT = 2
} Node_status;

#endif

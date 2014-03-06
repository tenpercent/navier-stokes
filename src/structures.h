typedef struct {
  double time_upper_boundary;     // time boundaries: [0, T]
  double space_upper_boundary;    // space boundaries: [0, X]
  double p_ro;                    // p(ro) is the function from the problem statement
                                  // we consider it is linear
                                  // i.e. p(ro) = p_ro * ro
  double viscosity;               // viscosity; also look into the problem statement. MUUUUU
} Gas_parameters;

typedef struct {
  unsigned X_nodes;
  unsigned T_nodes;
  unsigned total_nodes;
  double X_step;
  double T_step;
} Grid;

typedef enum {
  LEFT = 0,
  MIDDLE = 1,
  RIGHT = 2
} Node_status;

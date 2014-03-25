typedef struct {
  /* public stuff */
  unsigned size;
  unsigned *indices;
  double *elements;
  /* for private use */
  unsigned filled_row;
  unsigned filled_element;
} Sparse_matrix;

void Sparse_matrix_Construct (
    Sparse_matrix * this,
    unsigned size,
    unsigned nzcount);

void Sparse_matrix_Append_element (
    Sparse_matrix * this,
    unsigned row,
    unsigned col,
    double value);

void Sparse_matrix_Finish_filling (
    Sparse_matrix * this);

void Sparse_matrix_Destruct (
    Sparse_matrix * this);

void Sparse_matrix_Apply_to_vector (
    Sparse_matrix * this,
    double * vector,
    double * newVector);

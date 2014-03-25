typedef struct {
  unsigned size;
  unsigned *indices;
  double *elements;
} Sparse_matrix;

void Sparse_matrix_Construct (
    Sparse_matrix * this,
    unsigned size,
    unsigned nzcount);

void Sparse_matrix_Destruct (
    Sparse_matrix * this);

void Sparse_matrix_Apply_to_vector (
    Sparse_matrix * this,
    double * vector,
    double * newVector);

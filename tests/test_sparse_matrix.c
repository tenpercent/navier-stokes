#include <assert.h>
#include <math.h>
#include "sparse_matrix.h"

#define D_EQUAL(d1, d2) (fabs((d1) - (d2)) < 1e-10)
#define ASSERT_EQUAL(d1, d2) assert(D_EQUAL(d1, d2))

void test_filling_and_apply(Sparse_matrix *matrix) {
    double vector[4], new_vector[4];

    new_vector[0] = 99.;
    new_vector[1] = 99.;
    new_vector[2] = 99.;
    new_vector[3] = 99.;

    Sparse_matrix_Construct(matrix, 4, 8);

    /* To make Apply_to_vector work */
    matrix->elements[0] = 0.;
    matrix->elements[1] = 0.;
    matrix->elements[2] = 0.;
    matrix->elements[3] = 0.;

    Sparse_matrix_Append_element(matrix, 0, 3, 1.);
    Sparse_matrix_Append_element(matrix, 0, 0, 3.);
    Sparse_matrix_Append_element(matrix, 1, 1, 2.);
    Sparse_matrix_Append_element(matrix, 1, 3, 1.);
    Sparse_matrix_Append_element(matrix, 2, 1, 3.);
    Sparse_matrix_Append_element(matrix, 2, 0, 1.);
    Sparse_matrix_Append_element(matrix, 3, 2, 2.);
    Sparse_matrix_Append_element(matrix, 3, 3, 1.);
    Sparse_matrix_Finish_filling(matrix);

    vector[0] = 1;
    vector[1] = 0;
    vector[2] = 3;
    vector[3] = 2;

    /* / 3 0 0 1 \/ 1 \   / 5 \
       | 0 2 0 1 || 0 | = | 2 |
       | 1 3 0 0 || 3 |   | 1 |
       \ 0 0 2 1 /\ 2 /   \ 8 / */

    Sparse_matrix_Apply_to_vector(matrix, vector, new_vector);

    ASSERT_EQUAL(new_vector[0], 5.);
    ASSERT_EQUAL(new_vector[1], 2.);
    ASSERT_EQUAL(new_vector[2], 1.);
    ASSERT_EQUAL(new_vector[3], 8.);

    Sparse_matrix_Clear(matrix);

    /* To make Apply_to_vector work */
    matrix->elements[0] = 0.;
    matrix->elements[1] = 0.;
    matrix->elements[2] = 0.;
    matrix->elements[3] = 0.;

    Sparse_matrix_Append_element(matrix, 0, 0, 1.);
    Sparse_matrix_Append_element(matrix, 1, 0, 2.);
    Sparse_matrix_Append_element(matrix, 2, 0, 3.);
    Sparse_matrix_Append_element(matrix, 2, 2, 2.);
    Sparse_matrix_Append_element(matrix, 3, 0, 4.);
    Sparse_matrix_Finish_filling(matrix);

    Sparse_matrix_Apply_to_vector(matrix, vector, new_vector);

    /* / 1 0 0 0 \/ 1 \   / 1 \
       | 2 0 0 0 || 0 | = | 2 |
       | 3 0 2 0 || 3 |   | 9 |
       \ 4 0 0 0 /\ 2 /   \ 4 / */

    ASSERT_EQUAL(new_vector[0], 1.);
    ASSERT_EQUAL(new_vector[1], 2.);
    ASSERT_EQUAL(new_vector[2], 9.);
    ASSERT_EQUAL(new_vector[3], 4.);

    Sparse_matrix_Destruct(matrix);
}

int main() {
    Sparse_matrix matrix;
    test_filling_and_apply(&matrix);
    return 0;
}

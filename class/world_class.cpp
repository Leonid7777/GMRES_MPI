#include <algorithm>
#include "world_class.h"

void 
Matvec::mat_vec(double* vec, double* res)
{
    static constexpr int block_size = 128;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < size_of_matrix; ++i) {
        for (int j = 0; j < width_of_matrix; j += block_size) {
            for (int jj = j; jj < std::min(j + block_size, width_of_matrix); ++jj) {
                res[i] += matrix[i * width_of_matrix + jj] * vec[i];
            }
        }
    }
}

double*
Matvec::get_matrix()
{
    return matrix;
}

int
Matvec::get_size()
{
    return size_of_matrix;
}

int 
Matvec::get_width()
{
    return width_of_matrix;
}

Matvec::Matvec(double* mtx, const int size, const int width)
{
    matrix = mtx;
    size_of_matrix = size;
    width_of_matrix = width;
}
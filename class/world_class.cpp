#include <algorithm>
#include <omp.h>
#include "world_class.h"

void 
Matvec::mat_vec(double* vec, double* res) const
{
    static constexpr int block_size = 128;

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < width_of_matrix; j += block_size) {
        int top = std::min(j + block_size, width_of_matrix);
        for(int i = 0; i < size_of_matrix; i++) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] += matrix[i * size_of_matrix + jj] * vec[i];
            }
        }
    }
}

double*
Matvec::get_matrix() const
{
    return matrix;
}

int
Matvec::get_size() const
{
    return size_of_matrix;
}

int 
Matvec::get_width() const
{
    return width_of_matrix;
}

Matvec::Matvec(double* mtx, const int size, const int width)
{
    matrix = mtx;
    size_of_matrix = size;
    width_of_matrix = width;
}
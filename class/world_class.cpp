#include <algorithm>


class Matvec
{
private:
    int size_of_matrix;
    int width_of_matrix;
    double* matrix;

public:
    Matvec(double* mtx, const int size, const int width)
    {
        matrix = mtx;
        size_of_matrix = size;
        width_of_matrix = width;
    }

    double*
    get_matrix()
    {
        return matrix;
    }

    int
    get_size()
    {
        return size_of_matrix;
    }

    int
    get_width()
    {
        return width_of_matrix;
    }

    static void
    mat_vec(double* mat, double* vec, double* res, int n, int m) {
        int block_size = 16;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < n; i += block_size) {
            for (int j = 0; j < m; j += block_size) {
                for (int ii = i; ii < std::min(i + block_size, n); ++ii) {
                    for (int jj = j; jj < std::min(j + block_size, m); ++jj) {
                        res[ii] += mat[ii * m + jj] * vec[jj];
                    }
                }
            }
        }
    }

};
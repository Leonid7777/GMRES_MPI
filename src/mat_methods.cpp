#include <vector>
#include <algorithm>
#include <omp.h>


void 
eye_mat_place(std::vector<double>& mat, int n, int m, int place) 
{
    for(int i = place; i < m; i++) {
        mat[i * n + i] = 1;
    }
}

void 
both_position_mat_mul(std::vector<double>& res, std::vector<double>& vec_res, int n, int m, int position, double c, double s)
{
    double val_1;

    #pragma omp for schedule(static)
    for(int j = 0; j < m; j++) {
        val_1 = c * res[j * n + position] + s * res[j * n + position + 1];
        res[j * n + position + 1] = - s * res[j * n + position] + c * res[j * n + position + 1];
        res[j * n + position] = val_1;

        val_1 = c * vec_res[j * n + position] + s * vec_res[j * n + position + 1];
        vec_res[j * n + position + 1] = - s * vec_res[j * n + position] + c * vec_res[j * n + position + 1];
        vec_res[j * n + position] = val_1;
    }

    val_1 = c * vec_res[m * n + position] + s * vec_res[m * n + position + 1];
    vec_res[m * n + position + 1] = - s * vec_res[m * n + position] + c * vec_res[m * n + position + 1];
    vec_res[m * n + position] = val_1;
}

void 
place_mat_vec(std::vector<double> mat, double* vec, int m, int size_of_matrix)
{
    static constexpr int block_size = 128;
    double* res = new double[m];

    for(int i = 0; i < m; i++) {
        res[i] = 0;
    }

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < m; j += block_size) {
        int top = std::min(j + block_size, m);
        for (int i = 0; i < m; ++i) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] += mat[i * size_of_matrix + jj] * vec[i];
            }
        }
    }

    for(int i = 0; i < m; i++) {
        vec[i] = res[i];
    }

    delete[] res;
}

void
vec_sub_mat(double* res, double* mat, double* vec, int m, int n)
{   
    n += 1;
    static constexpr int block_size = 128;
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < m; j += block_size) {
        int top = std::min(j + block_size, m);
        for (int i = 0; i < n; ++i) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] -= mat[i * m + jj] * vec[i];
            }
        }
    }
}
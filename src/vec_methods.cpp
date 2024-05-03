#include <cmath>
#include <omp.h>

void
scal_prod(double* res, double* vec, double* mat, int m, int n)
{
    #pragma omp parallel for schedule(static)
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j < m; j++) {
            res[i] += vec[j] * mat[i * m + j];
        }
    }
}

double
norm_vec(const double* vec, int n)
{
    double val = 0;
    #pragma omp parallel for reduction(+:val)
    for(int i = 0; i < n; i++) {
        val += vec[i] * vec[i]; 
    }

    return std::sqrt(val);
}

void
scal_vec(double* vec, int n, double val)
{
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
        vec[i] /= val;
    }
}

double
rotatation(double a, double b)
{
    return std::sqrt(a * a + b * b);
}
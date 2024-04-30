#include <cmath>


void 
mat_vec(double* mat, double* vec, double* res, int m, int n)
{
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res[j] += mat[i * m + j] * vec[i];
        }
    }
}

double
vec_to_vec(double* f_vec, double* s_vec, int n)
{
    double scal = 0;
    for(int i = 0; i < n; i++) {
        scal += f_vec[i] * s_vec[i];
    }
    return scal;
}

void
vec_sub_vec(double* f_vec, double* s_vec, int n, double val) 
{
    for(int i = 0; i < n; i++) {
        f_vec[i] -= val * s_vec[i];
    }
}

double
vec_norm(double* vec, int n)
{
    double val = 0;
    for(int i = 0; i < n; i++) {
        val += vec[i] * vec[i]; 
    }

    return std::sqrt(val);
}

void
vec_del(double* vec, double n, double val)
{
    for(int i = 0; i < n; i++) {
        vec[i] /= val;
    }
}

double
hypotenuse(double a, double b)
{
    return std::sqrt(a * a + b * b);
}
#include <vector>
#include <algorithm>


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
    double* res = new double[m];

    for(int i = 0; i < m; i++) {
        res[i] = 0;
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            res[j] += mat[size_of_matrix * i + j] * vec[i];
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
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j < m; j++) {
            res[j] -= vec[i] * mat[i * m + j];
        }
    }
}
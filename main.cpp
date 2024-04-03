#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>


void
matrix_make(double* matrix, int size_of_matrix)
{
    std::srand(std::time(nullptr));
    for(int j = 0; j < size_of_matrix; j++) {
        for(int i = 0; i < size_of_matrix; i++) {
            matrix[i * size_of_matrix + j] = -1000.0 + 2000.0 * double(std::rand()) / RAND_MAX;
            // matrix[i * size_of_matrix + j] = i + j * size_of_matrix + 1;
        }
    }
}

void
right_part_make(double* right_part, int size_of_matrix, double& norm)
{
    std::srand(std::time(nullptr));
    for(int i = 0; i < size_of_matrix; i++) {
        right_part[i] = -1000.0 + 2000.0 * double(std::rand()) / RAND_MAX;
        // right_part[i] = i + 1;
        norm += right_part[i] * right_part[i];
    }
    norm = std::sqrt(norm);
}

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

void 
eye_mat(double* mat, int n, int m) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) {
                mat[i * n + j] = 1;
            } else {
                mat[i * n + j] = 0;
            }
        }
    }
}

void
mat_mul(double* mat, double* res, int n, int m) {
    double* prom = new double[n * m];

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            prom[i + j * n] = 0;
            for(int p = 0; p < n; p++) {
                prom[i + j * n] += mat[i + p * n] * res[j * n + p];
            }
        }
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res[i + j * n] = prom[i + j * n];
        }
    }

    delete[] prom;
}

void
position_mat_mul(double* res, int n, int m, int position, double c, double s) {

    double val_1;

    for(int j = 0; j < m; j++) {
        val_1 = c * res[j * n + position] + s * res[j * n + position + 1];
        res[j * n + position + 1] = - s * res[j * n + position] + c * res[j * n + position + 1];
        res[j * n + position] = val_1;
    }
}

void
transpose_mat(double* mat, int n, int m)
{
    double* res = new double[m * n];

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            res[j * m + i] = mat[i * n + j];
        }
    }

    for(int i = 0; i < m * n; i++) {
        mat[i] = res[i];
    }
}

int 
main()
{
    int krylov_count = 100;
    int size_of_matrix = 100;

    double* right_part = new double[size_of_matrix];
    double* res = new double[size_of_matrix];
    double* matrix = new double[size_of_matrix * size_of_matrix];

    double* krylov_subspaces = new double[(krylov_count + 1) * size_of_matrix];

    double* H = new double[(krylov_count + 1) * krylov_count];
    double* Q = new double[(krylov_count + 1) * (krylov_count + 1)];
    double* R = new double[(krylov_count + 1) * (krylov_count + 1)];

    double norm = 0;

    eye_mat(Q, krylov_count + 1, krylov_count + 1);
    eye_mat(R, krylov_count + 1, krylov_count);

    matrix_make(matrix, size_of_matrix);
    right_part_make(right_part, size_of_matrix, norm);

    for(int i = 0; i < size_of_matrix; i++) {
        res[i] = 0;
        krylov_subspaces[i] = right_part[i] / norm;
    }

    for(int j = 0; j < krylov_count; j++) {

        // It's OK
        mat_vec(matrix, krylov_subspaces + j * size_of_matrix, krylov_subspaces + (j + 1) * size_of_matrix, size_of_matrix, size_of_matrix);

        // It's OK
        for(int i = 0; i <= j; i++) {

            H[j * (krylov_count + 1) + i] = vec_to_vec(krylov_subspaces + (j + 1) * size_of_matrix, krylov_subspaces + i * size_of_matrix, size_of_matrix);
        }

        // It's OK
        for(int i = 0; i <= j; i++) {
            vec_sub_vec(krylov_subspaces + (j + 1) * size_of_matrix, krylov_subspaces + i * size_of_matrix, size_of_matrix, H[j * (krylov_count + 1) + i]);
        }

        // It's OK
        H[j * (krylov_count + 1) + j + 1] = vec_norm(krylov_subspaces + (j + 1) * size_of_matrix, size_of_matrix);


        if(H[j * (krylov_count + 1) + j + 1] == 0) {
            break;
        } else {
            vec_del(krylov_subspaces + (j + 1) * size_of_matrix, size_of_matrix, H[j * (krylov_count + 1) + j + 1]);
        }

    }
    
    double denom, c, s;

    for(int i = 0; i < krylov_count; i++) {

        denom = std::sqrt(H[i * (krylov_count + 1) + i] * H[i * (krylov_count + 1) + i] + H[i * (krylov_count + 1) + i + 1] * H[i * (krylov_count + 1) + i + 1]);

        c =  H[i * (krylov_count + 1) + i] / denom;
        s =  H[i * (krylov_count + 1) + i + 1] / denom;

        position_mat_mul(H, krylov_count + 1, krylov_count, i, c, s);
        position_mat_mul(Q, krylov_count + 1, krylov_count + 1, i, c, s);
    }

    double* r_p = new double[krylov_count];
    double* y = new double[krylov_count];

    for(int i = 0; i < krylov_count; i++) {
        r_p[i] = Q[i] * norm;
    }

    for(int i = krylov_count - 1; i >= 0; i--) {
        y[i] = r_p[i];
        for(int j = i + 1; j < krylov_count; j++) {
            y[i] -= H[j * (krylov_count + 1) + i] * y[j];
        }
        y[i] /=  H[i * (krylov_count + 1) + i];
    }

    for(int i = 0; i < krylov_count; i++) {
        for(int j = 0; j < size_of_matrix; j++) {
            res[j] += y[i] * krylov_subspaces[i * size_of_matrix + j];
        }
    }

    double* h_i_g = new double[size_of_matrix];
    mat_vec(matrix, res, h_i_g, size_of_matrix, size_of_matrix);

    double normis = 0;
    double n_r_p = 0;

    for(int j = 0; j < size_of_matrix; j++) {
        normis += (h_i_g[j] - right_part[j]) * (h_i_g[j] - right_part[j]);
        n_r_p += right_part[j] * right_part[j];
    }

    std::cout << std::sqrt(normis / n_r_p) << std::endl;
    

    return 0;
}
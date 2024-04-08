#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <vector>


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
eye_mat(std::vector<double>& mat, int n, int m) {
    // for(int i = 0; i < m; i++) {
    //     for(int j = 0; j < n; j++) {
    //         if(i == j) {
    //             mat[i * n + j] = 1;
    //         } else {
    //             mat[i * n + j] = 0;
    //         }
    //     }
    // }

    for(int i = 0; i < m; i++) {
        mat[i * n + i] = 1;
    }

}

void 
eye_mat_vec_place(std::vector<double>& mat, int n, int m, int place) {

    // for(int i = place; i < m; i++) {

    //     for(int j = 0; j < n; j++) {
    //         if(i == j) {
    //             mat[i * n + j] = 1;
    //         } else {
    //             mat[i * n + j] = 0;
    //         }
    //     }
    // }

    for(int i = place; i < m; i++) {
        mat[i * n + i] = 1;
    }

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
position_mat_mul_vec(std::vector<double>& res, int n, int m, int position, double c, double s) {
    double val_1;

    for(int j = 0; j < m; j++) {
        val_1 = c * res[j * n + position] + s * res[j * n + position + 1];
        res[j * n + position + 1] = - s * res[j * n + position] + c * res[j * n + position + 1];
        res[j * n + position] = val_1;
    }
}

void 
both_position_mat_mul(double* res, std::vector<double>& vec_res, int n, int m, int position, double c, double s)
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
        for(int j = 0; j < m; j++) {
            res[i] += mat[i + size_of_matrix * j] * vec[j];
        }
    }

    for(int i = 0; i < m; i++) {
        vec[i] = res[i];
    }

    delete[] res;
}

int 
main()
{
    int krylov_count = 0;
    int size_of_matrix = 100;
    double norm = 0;

    double* matrix = new double[size_of_matrix * size_of_matrix];
    double* right_part = new double[size_of_matrix];
    double* res = new double[size_of_matrix];

    matrix_make(matrix, size_of_matrix);
    right_part_make(right_part, size_of_matrix, norm);

    double* krylov_subspaces = new double[(size_of_matrix + 1) * size_of_matrix];
    double* H = new double[(size_of_matrix + 1) * size_of_matrix];

    int kr_v = std::max(2, (size_of_matrix + 1) / 10);
    std::vector<double> Q_vec ((size_of_matrix + 1) * kr_v);

    eye_mat(Q_vec, size_of_matrix + 1, kr_v);

    for(int i = 0; i < size_of_matrix; i++) {
        res[i] = 0;
        krylov_subspaces[i] = right_part[i] / norm;
    }

    double denom, c, s;
    double err = 10;

    while(err >= 0.001 && krylov_count < size_of_matrix) {

        if(krylov_count + 1 == kr_v) {
            int place = kr_v;
            kr_v = std::min(2 * kr_v, size_of_matrix + 1);
            Q_vec.resize(kr_v * (size_of_matrix + 1));
            eye_mat_vec_place(Q_vec, size_of_matrix + 1, kr_v, place);
        }

        mat_vec(matrix, krylov_subspaces + krylov_count * size_of_matrix, krylov_subspaces + (krylov_count + 1) * size_of_matrix, size_of_matrix, size_of_matrix);

        for(int i = 0; i <= krylov_count; i++) {
            H[krylov_count * (size_of_matrix + 1) + i] = vec_to_vec(krylov_subspaces + (krylov_count + 1) * size_of_matrix, krylov_subspaces + i * size_of_matrix, size_of_matrix);
        }

        for(int i = 0; i <= krylov_count; i++) {
            vec_sub_vec(krylov_subspaces + (krylov_count + 1) * size_of_matrix, krylov_subspaces + i * size_of_matrix, size_of_matrix, H[krylov_count * (size_of_matrix + 1) + i]);
        }

        H[krylov_count * (size_of_matrix + 1) + krylov_count + 1] = vec_norm(krylov_subspaces + (krylov_count + 1) * size_of_matrix, size_of_matrix);

        if(H[krylov_count * (size_of_matrix + 1) + krylov_count + 1] == 0) {
            krylov_count++;
            break;
        } else {
            vec_del(krylov_subspaces + (krylov_count + 1) * size_of_matrix, size_of_matrix, H[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);
        }

        place_mat_vec(Q_vec, H + krylov_count * (size_of_matrix + 1), krylov_count + 1, size_of_matrix + 1);

        denom = std::sqrt(H[krylov_count * (size_of_matrix + 1) + krylov_count] * H[krylov_count * (size_of_matrix + 1) + krylov_count] + H[krylov_count * (size_of_matrix + 1) + krylov_count + 1] * H[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);

        c =  H[krylov_count * (size_of_matrix + 1) + krylov_count] / denom;
        s =  H[krylov_count * (size_of_matrix + 1) + krylov_count + 1] / denom;

        // position_mat_mul(H, size_of_matrix + 1, krylov_count + 1, krylov_count, c, s);
        // position_mat_mul_vec(Q_vec, size_of_matrix + 1, krylov_count + 2, krylov_count, c, s);

        both_position_mat_mul(H, Q_vec, size_of_matrix + 1, krylov_count + 1, krylov_count, c, s);

        err = std::abs(Q_vec[krylov_count + 1]);
        krylov_count++;
    }

    std::cout << krylov_count << std::endl;

    double* y = new double[size_of_matrix];

    for(int i = krylov_count - 1; i >= 0; i--) {
        y[i] = Q_vec[i] * norm;
        for(int j = i + 1; j < krylov_count; j++) {
            y[i] -= H[j * (size_of_matrix + 1) + i] * y[j];
        }
        y[i] /=  H[i * (size_of_matrix + 1) + i];
    }

    for(int i = 0; i < krylov_count; i++) {
        for(int j = 0; j < size_of_matrix; j++) {
            res[j] += y[i] * krylov_subspaces[i * size_of_matrix + j];
        }
    }

    double* h_i_g = new double[size_of_matrix];
    mat_vec(matrix, res, h_i_g, size_of_matrix, size_of_matrix);

    double normis = 0;

    for(int j = 0; j < size_of_matrix; j++) {
        normis += (h_i_g[j] - right_part[j]) * (h_i_g[j] - right_part[j]);
    }

    std::cout << std::sqrt(normis) / norm << std::endl;

    delete[] right_part;
    delete[] res;
    delete[] matrix;
    delete[] krylov_subspaces;
    delete[] H;
    delete[] y;
    delete[] h_i_g;

    return 0;
}
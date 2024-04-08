#include <vector>
#include <cmath>
#include "vec_methods.h"
#include "mat_methods.h"


void
GMRES(int size_of_matrix, double* matrix, double* right_part, double* res)
{
    double norm = vec_norm(right_part, size_of_matrix);
    int krylov_count = 0;

    double* krylov_subspaces = new double[(size_of_matrix + 1) * size_of_matrix];
    double* H = new double[(size_of_matrix + 1) * size_of_matrix];

    int kr_v = std::max(2, (size_of_matrix + 1) / 10);
    std::vector<double> Q_vec ((size_of_matrix + 1) * kr_v);

    eye_mat_place(Q_vec, size_of_matrix + 1, kr_v, 0);

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
            eye_mat_place(Q_vec, size_of_matrix + 1, kr_v, place);
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

        both_position_mat_mul(H, Q_vec, size_of_matrix + 1, krylov_count + 1, krylov_count, c, s);

        err = std::abs(Q_vec[krylov_count + 1]);
        krylov_count++;
    }

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

    delete[] krylov_subspaces;
    delete[] H;
    delete[] y;
}
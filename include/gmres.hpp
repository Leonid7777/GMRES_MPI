#include <vector>
#include <cmath>
#include "vec_methods.h"
#include "mat_methods.h"
#include "../class/world_class.h"


template <class T>
void
GMRES(T& A, const double* right_part, double* res)
{
    int size_of_matrix = A.get_size();
    double norm = vec_norm(right_part, size_of_matrix);
    int krylov_count = 0;

    double* krylov_subspaces = new double[(size_of_matrix + 1) * size_of_matrix];

    int kr_v = std::max(2, (size_of_matrix + 1) / 10);
    std::vector<double> Q_vec ((size_of_matrix + 1) * kr_v);
    std::vector<double> H_vec ((size_of_matrix + 1) * (kr_v - 1));

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
            H_vec.resize((kr_v - 1) * (size_of_matrix + 1));
            eye_mat_place(Q_vec, size_of_matrix + 1, kr_v, place);
        }

        A.mat_vec(krylov_subspaces + krylov_count * size_of_matrix, krylov_subspaces + (krylov_count + 1) * size_of_matrix);

        scal_prod(&H_vec[krylov_count * (size_of_matrix + 1)], krylov_subspaces + (krylov_count + 1) * size_of_matrix, krylov_subspaces, size_of_matrix, krylov_count);

        vec_sub_mat(krylov_subspaces + (krylov_count + 1) * size_of_matrix, krylov_subspaces, &H_vec[krylov_count * (size_of_matrix + 1)], size_of_matrix, krylov_count);

        H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] = vec_norm(krylov_subspaces + (krylov_count + 1) * size_of_matrix, size_of_matrix);

        if(H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] == 0) {
            krylov_count++;
            break;
        } else {
            vec_del(krylov_subspaces + (krylov_count + 1) * size_of_matrix, size_of_matrix, H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);
        }

        place_mat_vec(Q_vec, &H_vec[krylov_count * (size_of_matrix + 1)], krylov_count + 1, size_of_matrix + 1);

        denom = hypotenuse(H_vec[krylov_count * (size_of_matrix + 1) + krylov_count], H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);

        c =  H_vec[krylov_count * (size_of_matrix + 1) + krylov_count] / denom;
        s =  H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] / denom;

        both_position_mat_mul(H_vec, Q_vec, size_of_matrix + 1, krylov_count + 1, krylov_count, c, s);

        err = std::abs(Q_vec[krylov_count + 1]);
        krylov_count++;
    }

    double* y = new double[size_of_matrix];

    for(int i = krylov_count - 1; i >= 0; i--) {
        y[i] = Q_vec[i] * norm;
        for(int j = i + 1; j < krylov_count; j++) {
            y[i] -= H_vec[j * (size_of_matrix + 1) + i] * y[j];
        }
        y[i] /=  H_vec[i * (size_of_matrix + 1) + i];
    }

    for(int i = 0; i < krylov_count; i++) {
        for(int j = 0; j < size_of_matrix; j++) {
            res[j] += y[i] * krylov_subspaces[i * size_of_matrix + j];
        }
    }

    delete[] krylov_subspaces;
    delete[] y;
}
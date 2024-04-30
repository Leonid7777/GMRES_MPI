#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "create_mat_r.h"
#include "vec_methods.h"
#include "gmres.h"


int 
main()
{
    double norm = 0;
    int size_of_matrix = 10;

    double* matrix = new double[size_of_matrix * size_of_matrix];
    double* right_part = new double[size_of_matrix];
    double* res = new double[size_of_matrix];

    matrix_make(matrix, size_of_matrix);
    right_part_make(right_part, size_of_matrix, norm);

    auto start = std::chrono::high_resolution_clock::now();

    GMRES(size_of_matrix, matrix, right_part, res);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Время работы программы: " << duration.count() << " миллисекунд" << std::endl;

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
    delete[] h_i_g;

    return 0;
}
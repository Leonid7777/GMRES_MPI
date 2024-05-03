#include <algorithm>


class Matvec
{
private:
    int size_of_matrix;
    int width_of_matrix;
    double* matrix;

public:
    Matvec(double* mtx, const int size, const int width);

    double* get_matrix() const;

    int get_size() const;

    int get_width() const;

    void mat_vec(double* vec, double* res) const;

};
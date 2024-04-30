void eye_mat_place(std::vector<double>& mat, int n, int m, int place);
void place_mat_vec(std::vector<double> mat, double* vec, int m, int size_of_matrix);
void both_position_mat_mul(std::vector<double>& res, std::vector<double>& vec_res, int n, int m, int position, double c, double s);
void vec_sub_mat(double* res, double* mat, double* vec, int m, int n);
#pragma once
#include "math/solver.hpp"
#include <vector>

#define MU_1 1
#define MU_2 1
#define ALPHA_SQ 0.1

#define U_X_0(x, y, t) cos(MU_2 * y) * exp(-(MU_1 * MU_1 + MU_2 * MU_2) * ALPHA_SQ * t)
#define U_X_1(x, y, t) -1 * cos(MU_2 * y) * exp(-(MU_1 * MU_1 + MU_2 * MU_2) * ALPHA_SQ * t)
#define U_Y_0(x, y, t) cos(MU_1 * x) * exp(-(MU_1 * MU_1 + MU_2 * MU_2) * ALPHA_SQ * t)
#define U_Y_1(x, y, t) -1 * cos(MU_1 * x) * exp(-(MU_1 * MU_1 + MU_2 * MU_2) * ALPHA_SQ * t)
#define U_XY(x, y) cos(MU_1 * x) * cos(MU_2 * y)

#define L 3.141593
#define T 1.0

#define N 6
#define M 6
#define K 12


#define H_X L / (N - 1)
#define H_Y L / (M - 1)
#define TAU T / K

//std::vector<Matrix> alterning_direction_method();
//Matrix step(Matrix prev);
//Matrix init();

void embed_vector(Matrix& A, Vector& b, bool cols, int ind);
Matrix horizontal(Matrix& prev_sol, float k);
Matrix vertical(Matrix& prev_sol, float k);
Matrix init();
Matrix step(Matrix& prev, int k);
std::vector<Matrix> alterning_direction_method();
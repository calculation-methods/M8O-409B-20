#pragma once
#include "math/solver.hpp"

#define U_X_0(x, y) y
#define U_X_1(x, y) 1 + y
#define U_Y_0(x, y) x
#define U_Y_1(x, y) 1 + x

#define L 1.0
#define N 11
#define M 11
#define H_X L / (N - 1)
#define H_Y L / (M - 1)

Matrix libman(float acc);
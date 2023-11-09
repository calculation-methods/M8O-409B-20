#pragma once

#include <fstream>
#include <cmath>
#include <vector>
#include "math/solver.hpp"

using namespace std;


#define ALPHA 0.5
#define LENGTH 3.141592653589
#define TIME_SEGMENT 2.0
#define N 101
#define M 101

const double h = LENGTH / (N - 1);
const double tau = TIME_SEGMENT / (M - 1);
const double c = tau * tau * ALPHA * ALPHA / (h * h);

double u_0(double x);
double du_0(double x);
double b_0(double t);
double b_1(double t);

void save_to_file(vector<Vector> sol);
Vector step_explicit(Vector x_0, Vector x_1, int k);
Vector step_implicit(Vector x_0, Vector x_1, int k);
Vector first_layer(Vector zero_layer);
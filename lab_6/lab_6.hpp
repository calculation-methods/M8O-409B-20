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

double u_0(double x);
double du_0(double x);
double b_0(double t);
double b_1(double t);

void save_to_file(vector<Vector> sol);
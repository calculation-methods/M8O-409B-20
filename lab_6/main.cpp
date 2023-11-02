#include <vector>
#include <iostream>
#include <cmath>
#include "solver.hpp"

#define ALPHA 0.8
#define LENGTH 3.1415
#define TIME_SEGMENT 1.0
#define N 11
#define M 11

using namespace std;

const double h = LENGTH / (N - 1);
const double tau = TIME_SEGMENT / (M - 1);

double u_0(double x)
{
    return sin(x);
}

double du_0(double x)
{
    return -ALPHA * sin(x);
}

double b_0(double t)
{
    return -sin(ALPHA * t);
}

double b_1(double t)
{
    return sin(ALPHA * t);
}

Vector first_layer(Vector zero_layer)
{
    Vector first_layer(N);

    Vector b(N);
    Matrix A(N, N);

    A[0][0] = 1;
    b[0] = b_0(tau);

    double c = tau * tau * ALPHA * ALPHA / (h * h);

    for(int i=1; i < N-1; i++)
    {
        
        b[i] = -1*(zero_layer[i] + du_0(h*i) * tau);
        A[i][i-1] = c / 2;
        A[i][i] = -(1 + c);
        A[i][i + 1] = c / 2;
    }
    A[N-1][N-1] = 1;
    b[N-1] = b_1(tau);

    return Tridiagonal_matrix_algorithm(A, b);
}
int main()
{   /*
    Matrix A(3, 3);
    A[0][0] = 2;
    A[0][1] = 1;

    A[1][0] = 0;
    A[1][1] = 2;
    A[1][2] = 0;

    A[2][1] = 0;
    A[2][2] = 2;

 
    Vector b(3);
    b[0] = 3;
    b[1] = 0;
    b[2] = 3;
    
    auto sol = Tridiagonal_matrix_algorithm(A, b);

    std::cout << sol[0] << " " <<sol[1] << " " <<sol[2] << std::endl;
    */
    //std::vector<Vector> sol;
    //Vector b(3);
    //b.get(1)= 8;
    //sol.push_back(b);
    //Vector c(3);
    //c.get(1)=8;
    //sol.push_back(c);
    //cout << sol[1][1];

    Vector zero_layer(N);
    for(int i = 0; i < N; i++)
        zero_layer[i] = 0.5;

    Vector first = first_layer(zero_layer);
    for(int i = 0; i < N; i++)
       cout << first[i] << endl;

    return 0;
}
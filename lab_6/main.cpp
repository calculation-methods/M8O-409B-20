#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "math/solver.hpp"

#define ALPHA 0.5
#define LENGTH 3.141592653589
#define TIME_SEGMENT 2.0
#define N 101
#define M 101

using namespace std;

const double h = LENGTH / (N - 1);
const double tau = TIME_SEGMENT / (M - 1);
const double c = tau * tau * ALPHA * ALPHA / (h * h);

double u_0(double x)
{
    return sin(x);
}

double du_0(double x)
{
    return -ALPHA * cos(x);
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
    /*
    Vector first_layer(N);

    for(int i=0; i<N; i++)
    {
        first_layer[i] = zero_layer[i] + tau * du_0(i*h);
    }
    */
    Vector b(N);
    Matrix A(N, N);

    A[0][0] = 1;
    b[0] = b_0(tau);

    //filling A & b
    
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

Vector step_explicit(Vector x_0, Vector x_1, int k)
{
    Vector x(N);

    x[0] = b_0(k*tau);

    for(int i=1; i < N-1; i++)
        x[i] = 2 * (1 - c) * x_1[i] + c * x_1[i-1] + c * x_1[i+1] - x_0[i];
    x[N-1] = b_1(k*tau);
    return x;
}

void save_to_file(vector<Vector> sol)
{
    std::ofstream file("sol.csv");
    for(int k=0; k<sol.size(); k++)
    {
        for(int j=0; j<sol[k].len; j++)
        {
            file << sol[k][j];
            if(j!=sol[k].len-1)
             file << ",";
        }
        file << "\n";
    }
    file.close();
    return;
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

    Vector zero(N);
    for(int i = 0; i < N; i++)
        zero[i] = u_0(i*h);

    Vector first = first_layer(zero);

    std::vector<Vector> sol;

    sol.push_back(zero);
    sol.push_back(first);

    for(int k=2; k<M; k++)
    {
        sol.push_back(step_explicit(sol[k-2], sol[k-1], k));
    }
    save_to_file(sol);
    return 0;
}
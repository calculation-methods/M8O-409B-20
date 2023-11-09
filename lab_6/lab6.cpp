#include "lab6.hpp"

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
    return -1 * sin(ALPHA * t);
}

double b_1(double t)
{
    return sin(ALPHA * t);
}

void save_to_file(vector<Vector> sol)
{
    ofstream file("sol.csv");
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

Vector step_implicit(Vector x_0, Vector x_1, int k)
{
    //Vector x(N);
    Matrix A(N, N);
    Vector b(N);
    //x[0] = b_0(k*tau);
    //x[N-1] = b_1(k*tau);
    A[0][0] = 1;
    b[0] = b_0(k*tau);//-2/(tau*tau) x_1[0] + 1/(tau*tau) X

    for(int i=1; i < N-1; i++)
    {
        A[i][i-1] = ALPHA*ALPHA / (h * h);
        A[i][i] = -(2*ALPHA*ALPHA / (h * h) + 1 / (tau*tau));
        A[i][i+1] = ALPHA*ALPHA / (h * h);
        b[i] = -2/(tau*tau) * x_1[i] + 1/(tau*tau) * x_0[i];
    }

    A[N-1][N-1] = 1;
    b[N-1] = b_1(k*tau);

    auto x = Tridiagonal_matrix_algorithm(A, b);
    return x;
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
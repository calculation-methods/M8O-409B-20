#include "matrix.hpp"
 
Matrix Tridiagonal_matrix_algorithm(Matrix A, Vector b)
{

    int n = A.rows;
    Vector alpha(n);
    Vector beta(n);
    Vector y(n);
    Vector x(n);
    
    // A[i] = A[i][i-1]
    // B[i] = A[i][i+1]
    // C[i] = A[i][i+1]
    y[0] = A[0][0];
    alpha[0] = -A[0][1] / y[0];
    beta[0] = b[0] / y[0];

    for(int i = 1; i < n-1; i++)
    {
        y[i] = alpha[i - 1] * A[i][i - 1] + A[i][i];
        alpha[i] = -A[i][i + 1] / y[i];
        beta[i] = (b[i] - A[i][i-1] * beta[i - 1]) / y[i];
    }
 
    y[n - 1] = alpha[n - 2] * A[n - 1][n - 2] + A[n - 1][n - 1];
    beta[n - 1] = (b[n - 1] - A[n-1][n-2] * beta[n - 2]) / y[n - 1];

    x[n-1] = beta[n-1];

    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}
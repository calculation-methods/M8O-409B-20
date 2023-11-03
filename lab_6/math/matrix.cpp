#include "matrix.hpp"


float& Matrix::get(int i, int j)
{
    return this->matr[i * cols + j];
}

float& Matrix::get(int i)
{
    return this->matr[i];
}

Matrix dot(Matrix A, Matrix B)
{
    if(A.cols != B.rows)
    {
        Matrix p(0, 0);
        return p;
    }
    
    Matrix P(A.rows, B.cols);
    for(int r = 0; r < A.rows; r++)
    {
        for(int c = 0; c < B.cols; c++)
        {
            float sum = 0;
            for(int ind = 0; ind < A.cols; ind++)
                sum += A.get(r, ind) * B.get(ind, c);
            P.get(r, c) = sum;
        }
    }
    return P;
}
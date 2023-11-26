#include "matrix.hpp"

//get matrix element
float const & Matrix::get(int i, int j) const
{
    return this->matr[i * cols + j];
}

//get vector element
float const & Vector::get(int i) const
{
    return this->matr[i];
}


//dot product = matrix multiplication
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
            P[r][c] = sum;
        }
    }
    return P;
}

//swap matrices fielda
void Matrix::swap(Matrix &other)
{
    std::swap(other.cols, cols);
    std::swap(other.rows, rows);
    std::swap(other.matr, matr);
}

//getting vect norm
float norm(Vector B)
{
    float max = -100000.0;
    for(int i=0; i<B.len; i++)
    {
        if(fabs(B.get(i)) > max)
        {
            max = fabs(B.get(i));
        }
    }
    return max;
}

void PrintMatrix(Matrix A)
{
    for(int i=0; i < A.rows; i++)
    {
        for(int j=0; j < A.cols; j++)
        {
            std::cout << A[i][j] << " ";
        }
        std:: cout << std::endl;
    }
    return;
}
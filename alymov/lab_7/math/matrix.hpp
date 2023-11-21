#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
//usual structures definitions - matrices and vectors

//matrices
struct Matrix
{

    void swap(Matrix &other);
    Matrix(int rows, int cols) : rows(rows), cols(cols), matr(new float[rows * cols]){}
    Matrix(Matrix const & M) : rows(M.rows), cols(M.cols), matr(new float[M.rows * M.cols])
    {
        for(int i = 0; i < M.rows; i++)
        {
            for(int j = 0; j < M.cols; j++)
            {   
                (*this)[i][j] = M.get(i, j);
            }
        }
    }


    float const & get(int i, int j) const;
    //float const & get(int i, int j);

    Matrix& operator=(Matrix other)
    {
        if(this != &other)
            Matrix(other).swap(*this);
        return *this;
    }

    Matrix operator+=(Matrix const & B)
    {
        if(B.rows != this->rows)
        {
            std::cout << "size mismatch! " << std::endl;
            return Matrix(B);
        }

        for(int i = 0; i < B.rows; i++)
        {
            for(int j = 0; j < B.cols; j++)
            {   
                (*this)[i][j] += B.get(i, j);
            }
        }
        return *this;
    }
  
    Matrix operator+(Matrix const & B)
    {
        return (*this)+=B;
    }

    Matrix operator*(float k)
    {
        Matrix mul(*this);
        for(int i = 0; i < this->rows; i++)
        {
            for(int j = 0; j < this->cols; j++)
            {   
                mul[i][j] *= k;
            }
        }
        return mul;
    }

    Matrix operator-=(Matrix const & B)
    {
        if(B.rows != this->rows)
        {
            std::cout << "size mismatch! " << std::endl;
            return Matrix(B);
        }
        for(int i = 0; i < B.rows; i++)
        {
            for(int j = 0; j < B.cols; j++)
            {   
                (*this)[i][j] -= B.get(i, j);
            }
        }
        return *this;
    }

    Matrix operator-(Matrix const & B)
    {
        Matrix copy(*this);
        return copy-=B;
    }

    float* operator[](int idx)
    {
        return ((this->matr) + idx * cols);
    }
    
    float* matr;
    unsigned rows;
    unsigned cols;

};

//vectors
struct Vector:Matrix
{
    Vector(int rows) : Matrix(rows, 1), len(rows){}

    Vector(Matrix A) : Matrix(A.rows, 1), len(A.rows){for(int i=0; i<A.rows; i++)*((this->matr) + i) = A[i][0];}

    float const & get(int i) const;

    float& operator[](int idx)
    {
        return *((this->matr) + idx);
    }
    using Matrix::operator=;
    
    //vector<double>
    unsigned len;
};

float norm(Vector B);
Matrix dot(Matrix A, Matrix B);
void PrintMatrix(Matrix A);
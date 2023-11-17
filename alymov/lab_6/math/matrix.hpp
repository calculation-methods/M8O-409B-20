#pragma once
struct Matrix
{

    Matrix(int rows, int cols) : rows(rows), cols(cols), matr(new float[rows * cols])
    {}

    float& get(int i, int j);
    float& get(int i);

    float* operator[](int idx)
    {
        return ((this->matr) + idx * cols);
    }

    //float& operator[](int )

    //float& operator[](int i);
    
    float *const matr;
    unsigned const rows;
    unsigned const cols;

};

struct Vector:Matrix
{
    Vector(int rows) : Matrix(rows, 1), len(rows)
    {}

    float& operator[](int idx)
    {
        return *((this->matr) + idx);
    }

    //vector<double>
    unsigned const len;
};

Matrix dot(Matrix A, Matrix B);
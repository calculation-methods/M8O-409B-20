#include <iostream>
#include "matrix.hpp"
#include "solver.hpp"

using namespace std;


int main()
{   
    Matrix A(3, 3);
    A[0][0] = 2;
    A[0][1] = 1;

    A[1][0] = 0;
    A[1][1] = 2;
    A[1][2] = 0;

    A[2][1] = 0;
    A[2][2] = 2;

 
    Vector b(3);
    b.get(0) = 3;
    b.get(1) = 3;
    b.get(2) = 3;
    
    auto sol = Tridiagonal_matrix_algorithm(A, b);

    std::cout << sol.get(0) << " " <<sol.get(1) << " " <<sol.get(2) << std::endl;
    return 0;
}
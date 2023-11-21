#include "math/solver.hpp"

int main()
{

    Vector b(2);
    //Vector c(3);
    //Matrix A(2, 2);

    b+=c;
    

    Vector sol = Jacobi(A, b);
    PrintMatrix(sol);
    return 0;
}
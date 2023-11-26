#include "lab_7.hpp"
#include <fstream>

//solution
void save_to_file(Matrix sol)
{
    std::ofstream file("sol.csv");
    for(int k=0; k<sol.rows; k++)
    {
        for(int j=0; j<sol.cols; j++)
        {
            file << sol[k][j];
            if(j!=sol.cols-1)
             file << ",";
        }
        file << "\n";
    }
    file.close();
    return;
}

int main(int argc, char* argv[])
{
    Matrix sol = libman(0.0001);
    PrintMatrix(sol);
    save_to_file(sol);
    return 0;
}
#include "lab_8.hpp"
#include <fstream>

void save_to_file(std::vector<Matrix>& system)
{
    std::ofstream file("sol.csv");

    for(int cnt=0; cnt<system.size(); cnt++){
        Matrix& sol = system[cnt];
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
    }
    file.close();
    return;
}
int main()
{
    std::vector<Matrix> sol = alterning_direction_method();
    PrintMatrix(sol[4]);
    save_to_file(sol);
    return 0;
}
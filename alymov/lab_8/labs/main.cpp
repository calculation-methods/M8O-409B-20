#include "lab_8.hpp"

int main(int arcg, char* argv[])
{
    if(argv[0] == "0")
        method = 0;
    else 
        method = 1;

    std::vector<Matrix> sol = alterning_direction_method();
    save_to_file(sol);
    return 0;
}
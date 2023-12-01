#include "lab_8.hpp"

int main()
{
    std::vector<Matrix> sol = alterning_direction_method();
    save_to_file(sol);
    return 0;
}
#include <iostream>
#include "lab6.hpp"

int main(int argc, char* argv[])
{  

    Vector zero(N);
    for(int i = 0; i < N; i++)
        zero[i] = u_0(i*h);

    Vector first = first_layer(zero);

    std::vector<Vector> sol;

    sol.push_back(zero);
    sol.push_back(first);

    if(argv[0] == "0")
    {
        for(int k=2; k<M; k++)
        {
            sol.push_back(step_implicit(sol[k-2], sol[k-1], k));
        }
    }
    else
    {
        for(int k=2; k<M; k++)
        {
            sol.push_back(step_explicit(sol[k-2], sol[k-1], k));
        }
    }

    save_to_file(sol);
    return 0;
}
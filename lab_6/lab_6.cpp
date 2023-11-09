#include "lab_6.hpp"

double u_0(double x)
{
    return sin(x);
}

double du_0(double x)
{
    return -ALPHA * cos(x);
}

double b_0(double t)
{
    return -1 * sin(ALPHA * t);
}

double b_1(double t)
{
    return sin(ALPHA * t);
}

void save_to_file(vector<Vector> sol)
{
    ofstream file("sol.csv");
    for(int k=0; k<sol.size(); k++)
    {
        for(int j=0; j<sol[k].len; j++)
        {
            file << sol[k][j];
            if(j!=sol[k].len-1)
             file << ",";
        }
        file << "\n";
    }
    file.close();
    return;
}
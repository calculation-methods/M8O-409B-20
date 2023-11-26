#include "lab_8.hpp"

std::vector<Matrix> alterning_direction_method(){
    std::vector<Matrix> sol;
    
    Matrix u(M, N);
    Matrix u_prev = init();

    sol.push_back(u_prev);

    for(int k = 0;k < K; k++){
        u = step(u_prev);
        sol.push_back(u);
        u = u_prev;
    }
    return sol;
}

Matrix init()
{
    Matrix initial(M, N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            initial[i][j] = U_XY(i * H_X, j * H_Y);
        }
    }
    return initial;
}

Matrix step(Matrix& prev)
{
    Matrix u(M, N);

    //for(int j)
}